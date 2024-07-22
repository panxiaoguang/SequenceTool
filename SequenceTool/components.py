import reflex as rx
from reflex.style import toggle_color_mode


def navbar() -> rx.Component:
    return rx.flex(
        rx.hstack(
            rx.hstack(
                rx.image(src="/DNA.jpeg", width="2.5em", height="auto"),
                rx.heading("SequenceTool", class_name="font-bold text-lg"),
                align_items="center"
            ),
            rx.hstack(
                rx.link(
                    "Home",
                    href="/",
                    color_scheme="gray"
                ),
                rx.link(
                    "AI Art",
                    href="/aiart",
                    color_scheme="gray"
                ),
                rx.link(
                    "AI Chat",
                    href="/aichat",
                    color_scheme="gray"
                ),
                spacing="5"
            ),
            align_items="center",
            spacing="3.5em",
            width="50%",
        ),
        rx.hstack(
            rx.button(
                rx.icon(tag="sun-medium"),
                on_click=toggle_color_mode,
                class_name=rx.color_mode_cond(
                    light="bg-white rounded-md text-slate-400", dark="bg-black text-white rounded-md")
            )
        ),
        class_name="border-b justify-between content-center p-3",
    )


def footer_item(text: str, href: str) -> rx.Component:
    return rx.link(rx.text(text, size="3"), href=href)


def footer_items_1() -> rx.Component:
    return rx.flex(
        rx.heading(
            "PRODUCTS", size="4", weight="bold", as_="h3"
        ),
        footer_item("Web Design", "/#"),
        footer_item("Web Development", "/#"),
        footer_item("E-commerce", "/#"),
        footer_item("Content Management", "/#"),
        footer_item("Mobile Apps", "/#"),
        spacing="4",
        text_align=["center", "center", "start"],
        flex_direction="column",
    )


def footer_items_2() -> rx.Component:
    return rx.flex(
        rx.heading(
            "RESOURCES", size="4", weight="bold", as_="h3"
        ),
        footer_item("Blog", "/#"),
        footer_item("Case Studies", "/#"),
        footer_item("Whitepapers", "/#"),
        footer_item("Webinars", "/#"),
        footer_item("E-books", "/#"),
        spacing="4",
        text_align=["center", "center", "start"],
        flex_direction="column",
    )


def social_link(icon: str, href: str) -> rx.Component:
    return rx.link(rx.icon(icon), href=href)


def socials() -> rx.Component:
    return rx.flex(
        social_link("instagram", "/#"),
        social_link("twitter", "/#"),
        social_link("facebook", "/#"),
        social_link("linkedin", "/#"),
        spacing="3",
        justify="end",
        width="100%",
    )


def footer() -> rx.Component:
    return rx.el.footer(
        rx.vstack(
            rx.flex(
                rx.vstack(
                    rx.hstack(
                        rx.image(
                            src="/DNA.jpeg",
                            width="3.25em",
                            height="auto",
                            border_radius="25%",
                        ),
                        rx.heading(
                            "SequenceTool",
                            size="7",
                            weight="bold",
                        ),
                        align_items="center",
                    ),
                    rx.text(
                        "© 2024 Xiaoguang, build with Reflex",
                        size="3",
                        white_space="nowrap",
                        weight="medium",
                    ),
                    spacing="4",
                    align_items=[
                        "center",
                        "center",
                        "start",
                    ],
                ),
                footer_items_1(),
                footer_items_2(),
                justify="between",
                spacing="6",
                flex_direction=["column", "column", "row"],
                width="100%",
            ),
            rx.divider(),
            rx.hstack(
                rx.hstack(
                    footer_item("Privacy Policy", "/#"),
                    footer_item("Terms of Service", "/#"),
                    spacing="4",
                    align="center",
                    width="100%",
                ),
                socials(),
                justify="between",
                width="100%",
            ),
            spacing="5",
            width="100%",
        ),
        width="100%",
        class_name=rx.color_mode_cond(light="p-4 bg-slate-100",
                                      dark="p-4 bg-neutral-800"),
    )
