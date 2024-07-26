import reflex as rx
from ..layout import template
from ..Markdown import mymarkdown, component_map
from ..utils import read_markdown_file


class BlogStat(rx.State):
    content: str

    async def get_blog(self):
        url = self.router.page.params.get("name", "no name")
        self.content = await read_markdown_file(f"{url}.md")
        yield


@rx.page(route="/Blog/[name]", on_load=BlogStat.get_blog)
@template
def blog() -> rx.Component:
    return rx.container(
        rx.flex(
            mymarkdown(BlogStat.content, component_map=component_map),
            width="80%",
            direction="column",
            class_name="mx-auto"
        )
    )
