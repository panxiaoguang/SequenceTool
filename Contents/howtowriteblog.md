# 如何用Reflex渲染Markdown文件生成博客

## Contents

## 实现思路

我们需要首先写一个渲染每一个blog页面的page函数，这个函数应该是一个动态路由，根据路由的URL，去Markdown文件所在路径寻找对应的文件并读取，然后将文件内容渲染为前端页面。

我们还需要一个能够列出所有blog的page函数，这个函数其实就是去Markdown路径，获取所有的文件名，并按照时间的顺序来排序。

## 首先把一些必须的功能性函数写出来

```python
## 计算时间戳
def get_file_create_time(filename):
    file = Path(os.path.join("Contents", filename))
    return datetime.datetime.fromtimestamp(file.stat().st_ctime).strftime("%Y-%m-%d")

## 读取Markdown的所有内容
async def read_markdown_file(filename):
    with open(os.path.join("Contents", filename), "r") as f:
        return f.read()

## 从内容里获取标题
def get_title(markdownContent):
    for line in markdownContent.splitlines():
        if line.strip().startswith("#"):
            return line.strip().replace("#", "").strip()

## 按照时间顺序列出所有Markdown文件
async def get_all_markdown_files():
    file_items = [[f, os.path.getctime(os.path.join(
        "Contents", f))] for f in os.listdir("Contents") if f.endswith(".md")]
    sorted_files = sorted(file_items, key=lambda x: x[1], reverse=True)
    return [file for file, _ in sorted_files]
```

## Blog 页面的写法

需要写一个on_load函数，即每次加载该页面，自动根据url获取对应的文件并读取渲染！！

```python
class BlogStat(rx.State):
    content: str

    async def get_blog(self):
        url = self.router.page.params.get("name", "no name")
        self.content = await read_markdown_file(f"{url}.md")
        yield


@rx.page(route="/Blog/[name]", on_load=BlogStat.get_blog)
def blog() -> rx.Component:
    return rx.container(
        rx.flex(
            mymarkdown(BlogStat.content, component_map=component_map),
            width="80%",
            direction="column",
            class_name="mx-auto"
        )
    )
```

## BlogView页面的写法

我们首先需要定义一个自定义变量，同时保存标题，时间和url路径；
然后同样需要一个on_load函数，去自动寻找所有的markdown文件并列出来。

```python
class BlogMeta(rx.Base):
    title: str
    url: str
    time: str


class BlogViewStat(rx.State):
    blogs: list[BlogMeta]

    async def get_all_blogs(self):
        self.blogs = []
        yield
        all_blogs = await get_all_markdown_files()
        for blog in all_blogs:
            fs = await read_markdown_file(blog)
            title = get_title(fs)
            url = f'/Blog/{blog.replace(".md", "")}'
            time = get_file_create_time(blog)
            self.blogs.append(BlogMeta(title=title, url=url, time=time))
```
前端页面就是用一个列表列出即可：

```python
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
def blogview() -> rx.Component:
    return rx.container(
        rx.box(
            rx.list.unordered(
                rx.foreach(BlogViewStat.blogs, list_blog),

            ),
            class_name='min-h-[50vh]'
        ),
    )
```

## 最后

这样，我们就可以实现自动读取blog的Markdown文件，并渲染出来了！