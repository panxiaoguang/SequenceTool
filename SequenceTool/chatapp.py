"""Welcome to Reflex! This file outlines the steps to create a basic app."""


import reflex as rx
from .layout import template
from openai import OpenAI


class QA(rx.Base):
    """A question and answer pair."""

    question: str
    answer: str


class ChatState(rx.State):
    """The app ChatState."""

    # A dict from the chat name to the list of questions and answers.
    chats: list[QA] = []

    # The current question.
    question: str

    # Whether we are processing the question.
    processing: bool = False

    models: list[str] = ["Qwen/Qwen2.5-7B-Instruct",
                         "Qwen/Qwen2.5-Coder-7B-Instruct",
                         "meta-llama/Meta-Llama-3.1-8B-Instruct"]

    # Default model
    model: str = "Qwen/Qwen2.5-7B-Instruct"

    apikey: str = rx.LocalStorage("")

    warning: bool = False

    def create_chat(self):
        """Create a new chat."""
        # Add the new chat to the list of chats.
        self.chats = []
        self.processing = False

    def close_warning(self):
        self.warning = False

    async def process_question(self, form_data: dict[str, str]):
        question = form_data["question"]
        if question == "":
            return
        if self.apikey != "":
            # Get the question from the form
            # Check if the question is empty

            model = self.openai_process_question

            async for value in model(question):
                yield value
        else:
            self.warning = True
            yield

    async def openai_process_question(self, question: str):
        """Get the response from the API.

        Args:
            form_data: A dict with the current question.
        """
        clients = OpenAI(
            base_url='https://api.siliconflow.cn/v1',
            api_key=self.apikey
        )
        # Add the question to the list of questions.
        qa = QA(question=question, answer="")
        self.chats.append(qa)

        # Clear the input and start the processing.
        self.processing = True
        yield

        # Build the messages.
        messages = [
            {
                "role": "system",
                "content": "You are a friendly chatbot named Reflex. Respond in markdown.",
            }
        ]
        for qa in self.chats:
            messages.append({"role": "user", "content": qa.question})
            messages.append({"role": "assistant", "content": qa.answer})

        # Remove the last mock answer.
        messages = messages[:-1]

        # Start a new session to answer the question.
        session = clients.chat.completions.create(
            model=self.model,
            messages=messages,
            stream=True,
        )

        # Stream the results, yielding after every word.
        for item in session:
            if hasattr(item.choices[0].delta, "content"):
                answer_text = item.choices[0].delta.content
                # Ensure answer_text is not None before concatenation
                if answer_text is not None:
                    self.chats[-1].answer += answer_text
                else:
                    # Handle the case where answer_text is None, perhaps log it or assign a default value
                    # For example, assigning an empty string if answer_text is None
                    answer_text = ""
                    self.chats[-1].answer += answer_text
                self.chats = self.chats
                yield

        # Toggle the processing flag.
        self.processing = False


message_style = dict(display="inline-block", padding="1em", border_radius="8px",
                     max_width=["30em", "30em", "50em", "50em", "50em", "50em"])


def set_apikey() -> rx.Component:
    return rx.dialog.root(
        rx.dialog.trigger(rx.button(rx.icon(tag="settings"),
                          variant="ghost", color_scheme="gray")),
        rx.dialog.content(
            rx.dialog.title("Set your API key"),
            rx.dialog.description(
                rx.text("Please find you API key from SiliconFlow:"),
                rx.link("https://cloud.siliconflow.cn/",
                        href="https://cloud.siliconflow.cn/"),
                size="2",
                margin_bottom="16px",
            ),
            rx.flex(
                rx.text(
                    "SiliconFlow API Key:",
                    as_="div",
                    size="2",
                    margin_bottom="4px",
                    weight="bold",
                ),
                rx.input(
                    value=ChatState.apikey,
                    placeholder="SiliconFlow API Key",
                    on_change=ChatState.set_apikey
                ),
                direction="column",
                spacing="3",
            ),
            rx.flex(
                rx.dialog.close(
                    rx.button(
                        "Cancel",
                        color_scheme="gray",
                        variant="soft",
                    ),
                ),
                rx.dialog.close(
                    rx.button("Save", on_click=ChatState.close_warning),
                ),
                spacing="3",
                margin_top="16px",
                justify="end",
            ),
        ),
    )


def message(qa: QA) -> rx.Component:
    """A single question/answer message.

    Args:
        qa: The question/answer pair.

    Returns:
        A component displaying the question/answer pair.
    """
    return rx.box(
        rx.box(
            rx.markdown(
                qa.question,
                background_color=rx.color("mauve", 4),
                color=rx.color("mauve", 12),
                **message_style,
            ),
            text_align="right",
            margin_top="1em",
        ),
        rx.box(
            rx.markdown(
                qa.answer,
                background_color=rx.color("accent", 4),
                color=rx.color("accent", 12),
                **message_style,
            ),
            text_align="left",
            padding_top="1em",
        ),
        width="100%",
    )


def chat() -> rx.Component:
    """List all the messages in a single conversation."""
    return rx.vstack(
        rx.box(rx.foreach(ChatState.chats, message), width="100%"),
        py="8",
        flex="1",
        width="100%",
        max_width="50em",
        padding_x="4px",
        align_self="center",
        overflow="hidden",
        padding_bottom="5em",
    )


def action_bar() -> rx.Component:
    """The action bar to send a new message."""
    return rx.center(
        rx.vstack(
            rx.hstack(
                rx.form(
                    rx.hstack(
                        rx.input(
                            rx.input.slot(
                                rx.tooltip(
                                    rx.icon("info", size=18),
                                    content="Enter a question to get a response.",
                                )
                            ),
                            placeholder="Type something...",
                            id="question",
                            width=["15em", "20em", "45em",
                                   "50em", "50em", "50em"],
                        ),
                        rx.button(
                            rx.cond(
                                ChatState.processing,
                                rx.spinner(size="3"),
                                rx.text("Send"),
                            ),
                            type="submit",
                            is_disabled=ChatState.processing,
                        ),
                        align_items="center",
                    ),
                    on_submit=ChatState.process_question,
                    reset_on_submit=True,
                ),
                rx.button(rx.text("New Chat"), on_click=ChatState.create_chat)
            ),
            rx.text(
                "ReflexGPT may return factually incorrect or misleading responses. Use discretion.",
                text_align="center",
                font_size=".75em",
                color=rx.color("mauve", 10),
            ),
            align_items="center",
        ),
        position="sticky",
        bottom="0",
        left="0",
        padding_y="16px",
        backdrop_filter="auto",
        backdrop_blur="lg",
        border_top=f"1px solid {rx.color('mauve', 3)}",
        background_color=rx.color("mauve", 2),
        align_items="stretch",
        width="100%",
    )


def select_model() -> rx.Component:
    return rx.center(
        rx.select(
            ChatState.models,
            default_value=ChatState.model,
            on_change=ChatState.set_model,
            variant="soft",
            radius="full",
            width="30%",
        ),
        set_apikey(),
        class_name="p-2",
        spacing="3"
    )


@rx.page("/aichat")
@template
def aichat() -> rx.Component:
    return rx.vstack(
        rx.cond(ChatState.warning,
                rx.callout("Access denied. Please set your API key!",
                           icon="triangle_alert",
                           color_scheme="red",
                           role="alert",),
                rx.flex(),
                ),
        select_model(),
        chat(),
        action_bar(),
        min_height="100vh",
        align_items="stretch",
        spacing="0",
    )
