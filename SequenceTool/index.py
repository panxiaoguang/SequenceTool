
import reflex as rx
from .utils import get_genome_seq, get_blat_data, parse_form_data, parse_blat_form, show_blat, rev_comp, to_lowercase, to_uppercase, global_alignment, render_codes
from Bio import SeqIO
import io
from .layout import template


class State(rx.State):
    sequence: str = ""
    seq_name: str = ""
    blat_data: list[list]
    modify_seq: str
    alignment_result: list[str]
    show_progress = False
    show_progress_blat = False

    async def getseq_submit(self, form_data: dict):
        form_data = parse_form_data(form_data)
        self.show_progress = True
        yield
        sequence_data = await get_genome_seq(**form_data)
        self.sequence = sequence_data['dna']
        self.seq_name = f">{sequence_data['chrom']}:{
            sequence_data['start']}-{sequence_data['end']}"
        self.show_progress = False
        yield

    async def blat_submit(self, form_data: dict):
        final_data = []
        form_data = parse_blat_form(form_data)
        self.show_progress_blat = True
        yield
        blat_data = await get_blat_data(**form_data)
        raw_blat_data = blat_data['blat']
        if form_data["Type"] == "protein":
            scale_factor = 3
        else:
            scale_factor = 1
        for raw_data in raw_blat_data:
            matches = int(raw_data[0])
            repMatches = int(raw_data[2])
            misMatches = int(raw_data[1])
            qinsert = int(raw_data[4])
            tinsert = int(raw_data[6])
            score = scale_factor*(matches+round(repMatches/2)) - \
                scale_factor*misMatches-qinsert-tinsert
            final_data.append([str(score), str(int(raw_data[11])+1), raw_data[12], raw_data[8],
                              raw_data[13], str(int(raw_data[15])+1), raw_data[16], raw_data[18]])
        self.blat_data = final_data
        self.show_progress_blat = False
        yield

    def revComplement(self):
        out_list = []
        if ">" in self.modify_seq:
            # input is fasta format and maybe more than one seqs
            for record in SeqIO.parse(io.StringIO(self.modify_seq), "fasta"):
                out_list.append(">"+str(record.id))
                out_list.append(str(record.seq.reverse_complement()))
            self.modify_seq = "\n".join(out_list)
        else:
            # no fasta format and only one seq
            self.modify_seq = rev_comp(self.modify_seq)

    def toUppercase(self):
        out_list = []
        if ">" in self.modify_seq:
            # input is fasta format and maybe more than one seqs
            for record in SeqIO.parse(io.StringIO(self.modify_seq), "fasta"):
                out_list.append(">"+str(record.id))
                out_list.append(str(record.seq).upper())
            self.modify_seq = "\n".join(out_list)
        else:
            self.modify_seq = to_uppercase(self.modify_seq)

    def toLowercase(self):
        out_list = []
        if ">" in self.modify_seq:
            # input is fasta format and maybe more than one seqs
            for record in SeqIO.parse(io.StringIO(self.modify_seq), "fasta"):
                out_list.append(">"+str(record.id))
                out_list.append(str(record.seq).lower())
            self.modify_seq = "\n".join(out_list)
        else:
            self.modify_seq = to_lowercase(self.modify_seq)

    def clear_modify_seq(self):
        self.modify_seq = ""

    def pairwise_alignment(self, form_data: dict):
        if form_data["end_gap"] == "false":
            aln, score = global_alignment(
                form_data["ref"], form_data["query"], int(form_data["gap_open"]), float(form_data["gap_extend"]), int(form_data["end_gap_open"]), float(form_data["end_gap_extend"]), False)
        else:
            aln, score = global_alignment(
                form_data["ref"], form_data["query"], int(form_data["gap_open"]), float(form_data["gap_extend"]), int(form_data["end_gap_open"]), float(form_data["end_gap_extend"]), True)
        self.alignment_result = aln.split("\n")
        self.alignment_result.append(f"Score: {score}")

    def remove_to_modify(self):
        self.modify_seq = self.sequence


def fetch_sequence() -> rx.Component:
    return rx.container(
        rx.card(
            rx.vstack(
                rx.flex(
                    rx.text("Fetch Sequence",
                            class_name="text-lg text-green-700"),
                    rx.tooltip(
                        rx.icon(tag="circle-dot", size=15),
                        content="This tool using UCSC API to fetch sequence from genome,input genome coordinates and select species to fetch sequence!",
                    ),
                    align="center"
                ),
                rx.form(
                    rx.vstack(
                        rx.hstack(
                            rx.text("Select Species:"),
                            rx.select(
                                ["Human",
                                 "Zebrafish",
                                 "Mouse"],
                                default_value="Human",
                                name="species"
                            ),
                            rx.text("Coordinates:"),
                            rx.input(placeholder="Genome coordinates",
                                     max_length=150,
                                     width="20em",
                                     name="coord"),
                        ),
                        rx.hstack(
                            rx.text("Select Direction:"),
                            rx.radio(
                                ["5'-3'", "3'-5'"],
                                direction="row",
                                spacing="2",
                                size="3",
                                name="reverse"
                            )
                        ),

                        rx.button("Fetch", type="submit"),

                    ),
                    on_submit=State.getseq_submit,
                    reset_on_submit=False,
                ),
                rx.button("send to modify",
                          on_click=State.remove_to_modify),
                rx.divider(),
                rx.box(
                    rx.cond(
                        State.show_progress,
                        rx.spinner(size="3"),
                        rx.flex(
                            rx.code(State.seq_name, variant="ghost"),
                            rx.code(State.sequence, variant="ghost"),
                            direction="column",
                        ),
                    ),
                    width="100%",
                    height="200px",
                    max_height="200px",
                    class_name=rx.color_mode_cond(
                        light="bg-slate-100 border rounded-md overflow-y-auto",
                        dark="bg-neutral-900 border rounded-md overflow-y-auto"
                    ),
                ),
                rx.divider(),
                rx.flex(
                    rx.text("Search Seq from Genome",
                            class_name="text-lg text-yellow-400"),
                    rx.tooltip(
                        rx.icon(tag="circle-dot", size=15),
                        content="This tool using UCSC blat API to search sequence from genome,input a short sequence to perform alignment, please don't use a sequence that is longer than 8Kb!",
                    ),
                    align="center"
                ),
                rx.form(
                    rx.vstack(
                        rx.hstack(
                            rx.text("Genome:"),
                            rx.select(
                                ["Human",
                                 "Zebrafish",
                                 "Mouse"],
                                default_value="Human",
                                name="genome"
                            ),
                            rx.text("Query type:"),
                            rx.select(
                                ["DNA",
                                 "protein",
                                 "translated RNA",
                                 "translated DNA"],
                                default_value="DNA",
                                name="query_type"
                            ),
                        ),
                        rx.text_area(
                            placeholder="Enter sequence to search",
                            name="seq",
                            rows='5',
                            width="100%"
                        ),
                        rx.button("Submit", type="submit")
                    ),
                    on_submit=State.blat_submit,
                ),
                rx.divider(),
                rx.box(
                    rx.cond(
                        State.show_progress_blat,
                        rx.spinner(size="3"),
                        rx.table.root(
                            rx.table.header(
                                rx.table.row(
                                    rx.table.column_header_cell("Score"),
                                    rx.table.column_header_cell("Start"),
                                    rx.table.column_header_cell("End"),
                                    rx.table.column_header_cell("Strand"),
                                    rx.table.column_header_cell("Chrom"),
                                    rx.table.column_header_cell("Start"),
                                    rx.table.column_header_cell("End"),
                                    rx.table.column_header_cell("Span"),
                                )
                            ),
                            rx.table.body(
                                rx.foreach(State.blat_data, show_blat)
                            )
                        )
                    ),
                    width="100%",
                    height="200px",
                    class_name=rx.color_mode_cond(
                        light="bg-slate-100 border rounded-md overflow-y-auto",
                        dark="bg-neutral-900 border rounded-md overflow-y-auto"
                    ),
                ),
                rx.divider(),
                rx.flex(
                    rx.text("Modify sequences",
                            class_name="text-lg text-violet-500"),
                    rx.tooltip(
                        rx.icon(tag="circle-dot", size=15),
                        content="This tool can reverse complement, to uppercase, to lowercase your sequence based on Biopython API!",
                    ),
                    align="center"
                ),
                rx.vstack(
                    rx.text_area(
                        value=State.modify_seq,
                        placeholder="Enter sequence to modify",
                        rows='5',
                        width="100%",
                        on_change=State.set_modify_seq,
                    ),
                    rx.hstack(
                        rx.button("Reverse Complement",
                                  on_click=State.revComplement),
                        rx.button("To Uppercase",
                                  on_click=State.toUppercase),
                        rx.button("To Lowercase",
                                  on_click=State.toLowercase),
                        rx.button("Clear",
                                  on_click=State.clear_modify_seq),
                    ),
                    width="100%"
                ),
                rx.divider(),
                rx.flex(
                    rx.text("DNA Pairwise Alignment",
                            class_name="text-lg text-pink-600"),
                    rx.tooltip(
                        rx.icon(tag="circle-dot", size=15),
                        content="This tool can perform global pairwise alignment based on Biopython API!",
                    ),
                    align="center"
                ),
                rx.form(
                    rx.vstack(
                        rx.text(
                            "Paste your sequence here, please don't add > before sequence"),
                        rx.text_area(
                            placeholder="Enter reference sequence",
                            name="ref",
                            rows='3',
                            width="100%"
                        ),
                        rx.text(
                            "Paste your sequence here, please don't add > before sequence"),
                        rx.text_area(
                            placeholder="Enter query sequence",
                            name="query",
                            rows='3',
                            width="100%"
                        ),
                        rx.hstack(
                            rx.text("End gap:"),
                            rx.select(
                                ['1', '5', '10', '15', '20', '25', '50', '100'],
                                default_value='10',
                                name="gap_open"
                            ),
                            rx.text("Gap extend:"),
                            rx.select(
                                ['0', '0.0005', '0.001', '0.05', '0.1', '0.2', '0.4',
                                    '0.5', '0.6', '0.8', '1.0', '5.0', '10.0'],
                                default_value='0.5',
                                name="gap_extend"
                            ),
                        ),
                        rx.flex(
                            rx.text("End gap:"),
                            rx.tooltip(
                                rx.icon("circle-dot", size=15),
                                content="false means not calculating end gap score!"
                            ),
                        ),
                        rx.select(
                            ["true", "false"],
                            default_value="false",
                            name="end_gap"
                        ),
                        rx.hstack(
                            rx.text("End gap open:"),
                            rx.select(
                                ['1', '5', '10', '15', '20', '25', '50', '100'],
                                default_value='10',
                                name="end_gap_open"
                            ),
                            rx.text("End gap extend:"),
                            rx.select(
                                ['0', '0.0005', '0.001', '0.05', '0.1', '0.2', '0.4',
                                    '0.5', '0.6', '0.8', '1.0', '5.0', '10.0'],
                                default_value='0.5',
                                name="end_gap_extend"
                            ),
                        ),
                        rx.button("Submit", type="submit")
                    ),
                    on_submit=State.pairwise_alignment,
                ),
                rx.divider(),
                rx.box(
                    rx.flex(
                        rx.foreach(State.alignment_result, render_codes),
                        direction="column"
                    ),
                    width="100%",
                    height="200px",
                    max_height="200px",
                    class_name=rx.color_mode_cond(
                        light="bg-slate-100 border rounded-md overflow-y-auto whitespace-pre",
                        dark="bg-neutral-900 border rounded-md overflow-y-auto whitespace-pre"
                    ),
                ),
            ),
        ),
    )


@rx.page("/")
@template
def index() -> rx.Component:
    return fetch_sequence()
