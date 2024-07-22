"""Welcome to Reflex! This file outlines the steps to create a basic app."""


import reflex as rx
from .layout import template
import aiohttp
import json

models_dict = {
    "QuWen-v1.5": "@cf/qwen/qwen1.5-14b-chat-awq",
    "OpenChat-v3.5": "@cf/openchat/openchat-3.5-0106",
    "LLAMA3": "@cf/meta/llama-3-8b-instruct"
}


class ChatState(rx.State):
    """The app state."""
    chat_history: list[tuple[str, str]]
    _messages = [{"role": "system", "content": "You are a smart and versatile artificial intelligence assistant. You will faithfully follow the user's instructions to solve and answer questions. Your answers are detailed and professional, and they do not contain offensive information."}]
    question: str
    model: str = "QuWen-v1.5"
    Generation: bool = False

    @rx.background
    async def get_answer(self):
        async with aiohttp.ClientSession() as session:
            answer = ""
            async with self:
                self._messages.append(
                    {"role": "user", "content": self.question})
                self.chat_history.append((self.question, answer))
                self.question = ""
            # Construct the post request to the API.
            model = models_dict[self.model]
            ### use your true account-id and api-key instead xxx
            async with session.post(
                f"https://api.cloudflare.com/client/v4/accounts/xxx/ai/run/{
                    model}",
                headers={
                    "Authorization": "Bearer xxx",
                    "Content-Type": "application/json"
                },
                json={
                    "messages": self.get_value(self._messages),
                    "stream": True,
                    "max_tokens": 2048
                }
            ) as response:
                if response.status == 200:
                    async for line in response.content:
                        async with self:
                            if not self.Generation:
                                return
                        decoded_line = line.decode('utf-8').strip()
                        if decoded_line.startswith('data: '):
                            json_data = decoded_line[len('data: '):]
                            # Handle the JSON data.
                            if decoded_line == 'data: [DONE]':
                                async with self:
                                    self.Generation = not self.Generation
                                    self._messages.append(
                                        {"role": "assistant", "content": answer})
                                    return
                            else:
                                data = json.loads(json_data)
                                answer += data.get('response')
                                async with self:
                                    self.chat_history[-1] = (
                                        self.chat_history[-1][0], answer)
                                    # push the final response to the role of system
                else:
                    async with self:
                        self.chat_history[-1] = (self.chat_history[-1][0],
                                                 f"Failed to fetch data: {response.status}")

    def clear_answer(self):
        self.chat_history = []
        self._messages = [{"role": "system", "content": "You are a smart and versatile artificial intelligence assistant. You will faithfully follow the user's instructions to solve and answer questions. Your answers are detailed and professional, and they do not contain offensive information."}]

    def toggle_running(self):
        self.Generation = not self.Generation
        if self.Generation:
            return ChatState.get_answer()


def qa(question: str, answer: str) -> rx.Component:
    component_style = {
        "h1": lambda text: rx.heading(
            text, size="5", margin_y="1em"
        ),
        "h2": lambda text: rx.heading(
            text, size="3", margin_y="1em"
        ),
        "h3": lambda text: rx.heading(
            text, size="1", margin_y="1em"
        ),
        "p": lambda text: rx.text(
            text, margin_y="1em"
        ),
        "code": lambda text: rx.code(text, color="orange", weight="bold"),
        "codeblock": lambda text, **props: rx.code_block(
            text, **props, theme="light", margin_y="1em", wrap_long_lines=True,
        ),
        "a": lambda text, **props: rx.link(
            text, **props, color="blue", _hover={"color": "orange"}
        ),
    }
    return rx.box(
        rx.box(
            rx.heading("You:", size="3", class_name="mt-1 font-semibold"),
            rx.markdown(question, component_map=component_style),
        ),
        rx.box(
            rx.heading("Bot:", size="3", class_name="mt-1 font-semibold"),
            rx.markdown(answer, component_map=component_style),
        ),
    )


def chat() -> rx.Component:
    return rx.container(
        rx.center(
            rx.select(
                [
                    "QuWen-v1.5",
                    "OpenChat-v3.5",
                    "LLAMA3"
                ],
                default_value=ChatState.model,
                radius="large",
                width="200px",
                on_change=ChatState.set_model,

            ),
            width="100%"
        ),
        rx.scroll_area(
            rx.center(
                rx.box(
                    rx.foreach(
                        ChatState.chat_history,
                        lambda messages: qa(messages[0], messages[1]),
                    ),
                    width='60%',
                ),
            ),
            class_name="mt-14",
            type='hover',
            height='64vh',
        ),
        rx.center(
            rx.text_area(
                placeholder="Ask a question",
                value=ChatState.question,
                on_change=ChatState.set_question,
                width="600px",
                auto_height=True,
            ),
            rx.button(
                rx.cond(
                    ~ChatState.Generation,
                    "Generate",
                    "Stop"
                ),
                height="4.5em",
                width="6.5em",
                variant="soft",
                radius="medium",
                _hover={"cursor": "pointer"},
                on_click=ChatState.toggle_running,
            ),
            rx.button(
                "New Chat",
                height="4.5em",
                width="6.5em",
                radius="medium",
                on_click=ChatState.clear_answer,
                _hover={"cursor": "pointer"},
                color_scheme="iris",
            ),
            spacing="1",
            width="100%",
            class_name="pt-4"
        )
    )


@rx.page("/aichat")
@template
def aichat() -> rx.Component:
    return chat()
