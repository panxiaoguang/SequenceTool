import reflex as rx
from .layout import template
from .utils import get_file_create_time, read_markdown_file, get_title, get_all_markdown_files


class BlogMeta(rx.Base):
    title: str
    url: str
    time: str


class BlogViewStat(rx.State):
    blogs: list[BlogMeta]
    show: bool = False

    async def get_all_blogs(self):
        self.blogs = []
        self.show = True
        yield
        all_blogs = await get_all_markdown_files()
        for blog in all_blogs:
            fs = await read_markdown_file(blog)
            title = get_title(fs)
            url = f'/Blog/{blog.replace(".md", "")}'
            time = get_file_create_time(blog)
            self.blogs.append(BlogMeta(title=title, url=url, time=time))
        self.show = False


def list_blog(Blog: BlogMeta) -> rx.Component:
    return rx.list_item(
        rx.flex(
            rx.text(Blog.time),
            rx.link(Blog.title, href=Blog.url),
            spacing='3',
            direction='row'
        )
    )


@rx.page("/Blog", on_load=BlogViewStat.get_all_blogs)
@template
def blogview() -> rx.Component:
    return rx.container(
        rx.box(
            rx.cond(
                BlogViewStat.show,
                rx.chakra.stack(
                    rx.chakra.skeleton_circle(size="30px"),
                    rx.chakra.skeleton_text(no_of_lines=10),
                    width="80%"
                ),
                rx.list.unordered(
                    rx.foreach(BlogViewStat.blogs, list_blog),

                ),
            ),
            class_name='min-h-[50vh]'
        ),
    )
