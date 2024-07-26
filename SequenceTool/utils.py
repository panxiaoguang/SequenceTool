import aiohttp
import reflex as rx
from Bio.Seq import reverse_complement
from Bio import Align
from Bio.Align import substitution_matrices
from pathlib import Path
import datetime
import os


async def fetch(session, url, params):
    async with session.get(url, params=params) as response:
        return await response.json()


async def fetch_blat(session, url, params):
    async with session.get(url, params=params) as response:
        return await response.json(content_type="text/html")


async def get_genome_seq(genome, chrom, start, end, reverse) -> dict:
    base_url = "https://api.genome.ucsc.edu/getData/sequence"
    input_params = {
        "genome": genome,
        "chrom": chrom,
        "start": start,
        "end": end,
        "revComp": reverse
    }
    async with aiohttp.ClientSession() as session:
        data = await fetch(session, base_url, input_params)
        return data


async def get_blat_data(userSeq, Type, db, output) -> dict:
    base_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    input_params = {
        "userSeq": userSeq,
        "type": Type,
        "db": db,
        "output": output
    }
    async with aiohttp.ClientSession() as session:
        data = await fetch_blat(session, base_url, input_params)
        return data


def parse_coord(coord_string):
    chrom, position = coord_string.split(":")
    start, end = position.split("-")
    return chrom, int(start), int(end)


def parse_form_data(formdata):
    newData = {}
    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData['genome'] = genome_dict[formdata["species"]]
    chrom, start, end = parse_coord(formdata["coord"])
    newData['chrom'] = chrom
    newData['start'] = start-1
    newData['end'] = end
    if formdata["reverse"] == "3'-5'":
        newData['reverse'] = 1
    else:
        newData['reverse'] = 0
    return newData


def parse_blat_form(formdata):
    newData = {}
    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData['db'] = genome_dict[formdata["genome"]]
    newData['Type'] = formdata['query_type']
    newData['userSeq'] = formdata['seq']
    newData['output'] = "json"
    return newData


def show_blat(blat_result_row: list):
    return rx.table.row(
        rx.table.cell(blat_result_row[0]),
        rx.table.cell(blat_result_row[1]),
        rx.table.cell(blat_result_row[2]),
        rx.table.cell(blat_result_row[3]),
        rx.table.cell(blat_result_row[4]),
        rx.table.cell(blat_result_row[5]),
        rx.table.cell(blat_result_row[6]),
        rx.table.cell(blat_result_row[7])
    )


def rev_comp(seq):
    return reverse_complement(seq)


def to_uppercase(seq):
    return seq.upper()


def to_lowercase(seq):
    return seq.lower()


def global_alignment(ref, query, gap_open, gap_extend, end_gap_open, end_gap_extend, end_gap):
    matrix = substitution_matrices.load("NUC.4.4")
    if not end_gap:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_gap_score=0
        )
    else:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_open_gap_score=0-end_gap_open,
            end_extend_gap_score=0-end_gap_extend
        )
    ref = ref.strip().upper()
    query = query.strip().upper()
    alignments = aligner.align(ref, query)
    return alignments[0].format(), alignments[0].score


def render_codes(code_row: str) -> rx.text:
    return rx.text(code_row, class_name="font-mono")

#################################################################################################################
############################ AI ART PAGE##########################################################################
#################################################################################################################


async def get_picture_CF(account_id, api_token, model, prompts, nprompts, w, h, nstep, cfg_scale):
    base_url = f"https://api.cloudflare.com/client/v4/accounts/{
        account_id}/ai/run/{model}"
    headers = {
        "Authorization": f"Bearer {api_token}",
        "Content-Type": "application/json"
    }
    input_params = {
        "guidance": cfg_scale,
        "height": h,
        "negative_prompt": nprompts,
        "num_steps": nstep,
        "width": w,
        "prompt": prompts
    }
    async with aiohttp.ClientSession() as session:
        async with session.post(base_url, json=input_params, headers=headers) as response:
            if response.status == 200:
                data = await response.read()
                return data
            else:
                return None


def get_file_create_time(filename):
    file = Path(os.path.join("Contents", filename))
    return datetime.datetime.fromtimestamp(file.stat().st_ctime).strftime("%Y-%m-%d")


async def read_markdown_file(filename):
    with open(os.path.join("Contents", filename), "r") as f:
        return f.read()


def get_title(markdownContent):
    for line in markdownContent.splitlines():
        if line.strip().startswith("#"):
            return line.strip().replace("#", "").strip()


async def get_all_markdown_files():
    file_items = [[f, os.path.getctime(os.path.join(
        "Contents", f))] for f in os.listdir("Contents") if f.endswith(".md")]
    sorted_files = sorted(file_items, key=lambda x: x[1], reverse=True)
    return [file for file, _ in sorted_files]
