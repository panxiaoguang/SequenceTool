import reflex as rx
from .utils import get_picture_SC
from .layout import template
from reflex_image_zoom import image_zoom
from openai import OpenAI


class Paint(rx.State):
    model: str = "FLUX.1-schnell"
    prompt: str = ""
    negtive_prompt: str = ""
    image_size: str = "1024x1024"
    steps: int = 20
    scale: float = 7.5
    translate: bool = False
    show_progress: bool = False
    show_advanced: bool = True
    image_apikey: str = rx.LocalStorage("")
    chat_endpoint: str = rx.LocalStorage("")
    chat_apikey: str = rx.LocalStorage("")
    chat_model: str = rx.LocalStorage("gpt-4o-mini")
    url: str = "/beauty.png"
    warning: bool = False

    async def paint(self):
        if self.image_apikey != "" and self.chat_apikey != "" and self.chat_endpoint != "":
            clients = OpenAI(
                base_url=self.chat_endpoint,
                api_key=self.chat_apikey
            )
            self.show_progress = True
            yield
            loading_picture, newpromt = await get_picture_SC(clients, self.model, self.image_apikey, self.prompt, self.negtive_prompt, self.image_size, self.steps, self.scale, self.translate, chatmodel=self.chat_model)
            if loading_picture is not None:
                self.url = loading_picture
                self.prompt = newpromt
            self.show_progress = False
            yield
        else:
            self.warning = True
            yield

    @rx.var
    def get_image_size_list(self) -> list[str]:
        # if self.model == "FLUX.1-schnell":
        #    return ["1024x1024", "512x1024", "768x512", "768x1024", "1024x576", "576x1024"]
        # else:
        #    return ["1024x1024", "1024x2048", "1536x1024", "1536x2048", "2048x1152", "1152x2048"]
        return ["1024x1024", "512x1024", "768x512", "768x1024", "1024x576", "576x1024",
                "1024x2048", "1536x1024", "1536x2048", "2048x1152", "1152x2048"]

    def change_translate(self, checked: bool):
        self.translate = checked

    def change_advanced(self):
        self.show_advanced = not self.show_advanced

    def set_steps(self, value: int):
        self.steps = value

    def set_scale(self, value: float):
        self.scale = value

    def close_warning(self):
        self.warning = False


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
                    value=Paint.image_apikey,
                    placeholder="SiliconFlow API Key",
                    on_change=Paint.set_image_apikey
                ),
                rx.text(
                    "Magic Prompt endpoint:",
                    as_="div",
                    size="2",
                    margin_bottom="4px",
                    weight="bold",
                ),
                rx.input(
                    value=Paint.chat_endpoint,
                    placeholder="https://api.openai.com/v1",
                    on_change=Paint.set_chat_endpoint
                ),
                rx.text(
                    "Magic Prompt API key:",
                    as_="div",
                    size="2",
                    margin_bottom="4px",
                    weight="bold",
                ),
                rx.input(
                    value=Paint.chat_apikey,
                    placeholder="xx-xxxxxxxxxxxxxxxxxxxxxxxxxx",
                    on_change=Paint.set_chat_apikey
                ),
                rx.text(
                    "Magic Prompt Model:",
                    as_="div",
                    size="2",
                    margin_bottom="4px",
                    weight="bold",
                ),
                rx.input(
                    value=Paint.chat_model,
                    placeholder="gpt-4o-mini",
                    on_change=Paint.set_chat_model
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
                    rx.button("Save", on_click=Paint.close_warning),
                ),
                spacing="3",
                margin_top="16px",
                justify="end",
            ),
        ),
    )


def paint_image():
    return rx.flex(
        rx.vstack(
            rx.cond(Paint.warning,
                    rx.callout("Access denied. Please set your API key!",
                               icon="triangle_alert",
                               color_scheme="red",
                               role="alert",),
                    rx.flex(),
                    ),
            rx.flex(
                rx.text("Model selection:"),
                set_apikey(),
                class_name="justify-between",
                width="100%",
            ),
            rx.select(
                ["FLUX.1-schnell", "FLUX.1-dev", "stable-diffusion-3-medium",
                    "stable-diffusion-xl-base-1.0"],
                default_value="FLUX.1-schnell",
                width="100%",
                on_change=Paint.set_model,
            ),
            rx.flex(
                rx.flex(rx.text("Prompt:"),),
                rx.flex(
                    rx.switch(checked=Paint.translate,
                              on_change=Paint.change_translate,
                              ),
                    rx.text("Maigc prompt"),
                    spacing="2",
                    align="center",
                ),
                class_name="justify-between",
                width="100%",
            ),
            rx.text_area(placeholder="What do you want to paint?",
                         size="3",
                         rows="5",
                         width="100%",
                         on_blur=Paint.set_prompt,
                         class_name="resize"
                         ),
            rx.text("Dimensions:"),
            rx.select(
                Paint.get_image_size_list,
                default_value="1024x1024",
                width="100%",
                on_change=Paint.set_image_size,
            ),
            rx.button("Advanced Options",
                      variant="ghost",
                      color_scheme="gray",
                      size="3",
                      on_click=Paint.change_advanced,),
            rx.cond(
                Paint.show_advanced,
                rx.flex(),
                rx.flex(
                    rx.text("Negative prompt:"),
                    rx.text_area(placeholder="What do you want to avoid?",
                                 size="3",
                                 rows="5",
                                 width="100%",
                                 on_blur=Paint.set_negtive_prompt,
                                 class_name="resize"
                                 ),
                    rx.flex(
                        rx.vstack(
                            rx.text("Steps:"),
                            rx.input(type="number",
                                     value=Paint.steps,
                                     on_change=Paint.set_steps,
                                     ),
                            width="40%",
                        ),
                        rx.vstack(
                            rx.text("Scale:"),
                            rx.input(type="number",
                                     value=Paint.scale,
                                     on_change=Paint.set_scale,
                                     ),
                            width="40%",
                        ),
                        width="100%",
                        spacing="2",
                    ),
                    width="100%",
                    direction="column",
                ),
            ),
            rx.button("Generate",
                      size="3",
                      width="100%",
                      color_scheme="iris",
                      radius="medium",
                      _hover={"cursor": "pointer"},
                      on_click=Paint.paint,
                      ),
            class_name="border-2 px-3 py-2",
            width="35%"),
        rx.flex(
            rx.cond(Paint.show_progress,
                    rx.spinner(size="3"),
                    rx.vstack(
                        image_zoom(rx.image(src=Paint.url,
                                            alt="Generated Image", width="50%")),
                        rx.divider(
                            class_name="border-dotted border-2 border-indigo-600"),
                        rx.text("Your Prmopt:", weight="bold"),
                        rx.text(
                            Paint.prompt, class_name="italic font-mono text-sm text-lime-500"),
                        width="100%", class_name="p-2",
                    ),),
            width="70%"))


@rx.page("/aiart")
@template
def aiart() -> rx.Component:
    return paint_image()
