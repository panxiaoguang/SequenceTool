import reflex as rx
from PIL import Image
from reflex_image_zoom import image_zoom
from .utils import get_picture_CF
from .layout import template
import io

model_dict = {
    "dreamshaper-8-lcm": "@cf/lykon/dreamshaper-8-lcm",
    "stable-diffusion-xl-base-1.0": "@cf/stabilityai/stable-diffusion-xl-base-1.0",
    "stable-diffusion-xl-lightning": "@cf/bytedance/stable-diffusion-xl-lightning"
}


class Paint(rx.State):
    model: str = "stable-diffusion-xl-base-1.0"
    prompt: str = "a cat"
    nprompt: str = "low quality"
    cfg_scale: float = 7.5
    nsteps: int = 20
    width: int = 1024
    height: int = 1024
    image = Image.open("default_picture/a_cat.png")
    show_progress: bool = False

    async def paint(self):
        self.show_progress = True
        yield
        loading_picture = await get_picture_CF("xxx", "xxx", model_dict[self.model], self.prompt, self.nprompt, self.width, self.height, self.nsteps, self.cfg_scale)
        if loading_picture is not None:
            self.image = Image.open(io.BytesIO(loading_picture))
        self.show_progress = False
        yield


def paint_image():
    return rx.container(
        rx.box(
            rx.heading("SD painting", size="5"),
            rx.flex(
                rx.text("Select a model"),
                rx.select(
                    ["dreamshaper-8-lcm",
                     "stable-diffusion-xl-base-1.0",
                     "stable-diffusion-xl-lightning"],
                    default_value=Paint.model,
                    on_change=Paint.set_model,
                ),
                spacing="2",
                align="center",
            ),
            rx.flex(
                rx.text("Positive prompt:"),
                rx.box(
                    rx.text_area(placeholder="eg: A photo of a cat",
                                 size="3",
                                 rows="5",
                                 on_blur=Paint.set_prompt),
                    width="100%",
                ),
                rx.text("Negative prompt:"),
                rx.box(
                    rx.text_area(placeholder="eg: ugly",
                                 size="3",
                                 rows="5",
                                 on_blur=Paint.set_nprompt),
                    width="100%",
                ),
                direction="column",
                spacing="2",
                width="100%",
            ),
            rx.flex(
                rx.text("CFG_sclae:"),
                rx.input(
                    type="number",
                    value=Paint.cfg_scale,
                    on_change=Paint.set_cfg_scale,
                ),
                rx.text("Steps:"),
                rx.input(
                    type="number",
                    value=Paint.nsteps,
                    on_change=Paint.set_nsteps,
                ),
                spacing="2",
                width="100%",
                align="center",
                class_name="mt-2"
            ),
            rx.flex(
                rx.text("width:"),
                rx.input(
                    type="number",
                    value=Paint.width,
                    on_change=Paint.set_width
                ),
                rx.text("height:"),
                rx.input(
                    type="number",
                    value=Paint.height,
                    on_change=Paint.set_height,
                ),
                rx.button(rx.icon(tag="dna"),
                          "Generate",
                          size="3",
                          color_scheme="iris",
                          radius="medium",
                          _hover={"cursor": "pointer"},
                          on_click=Paint.paint,
                          ),
                spacing="2",
                width="100%",
                align="center",
            ),
            rx.flex(rx.divider(size="2", width="100%"),
                    class_name="my-5",
                    width="100%",
                    ),
            rx.box(
                rx.text("Generated Image:"),
                rx.card(
                    rx.cond(
                        Paint.show_progress,
                        rx.spinner(size="3"),
                        image_zoom(rx.image(src=Paint.image,
                                            alt="Generated Image")),
                    ),
                    width="400px",
                    height="400px",
                    variant="classic",
                    class_name="mx-auto",
                ),
                width="100%",
                height="50%",

            ),
            class_name="mt-5 border-2 border-gray-300 p-5 rounded-md",
        ),
        class_name="mt-5",
    )


@rx.page("/aiart")
@template
def aiart() -> rx.Component:
    return paint_image()
