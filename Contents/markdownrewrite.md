# 如何为Reflex自带的markdown添加额外的插件

## Contents

Reflex是一个对python友好的网站app框架，本站的所有功能都是基于Reflex实现的。但是其自带的Markdown功能并不具备TOC(table of contents)的功能，所以为了让其生成的blog具有目录功能，并且点击目录可以跳转到指定位置，我们需要添加额外的插件。

## 查找需要的插件

通过查找[awesome-remark](https://github.com/remarkjs/awesome-remark)仓库，我们发现了[remark-toc](https://github.com/remarkjs/remark-toc)是我们需要的插件，但是它还需要配合另外的两个插件，分别是[rehype-slug](https://github.com/rehypejs/rehype-slug)和[rehype-autolink-headings](https://github.com/rehypejs/rehype-autolink-headings), 其中，`rehype-slug`的作用是为每个标题自动生成一个id, 而`rehype-autolink-headings`的作用这是同时为每个标题生成一个`<a/>`标签，从而实现点击TOC可以跳转对应标题的功能。

## 通过类函数重载的方法将插件导入

因为我们只是需要添加额外的插件，所以大部分功能无需改动，我们只需要修改原`Markdown`类的两个函数，分别是`_get_import`和`_render`。

```python
from reflex.utils.imports import ImportDict, ImportVar
from reflex.vars import Var
from reflex.components.markdown import Markdown
from reflex.components.tags.tag import Tag
from reflex.components.component import Component
import reflex as rx
_REMARK_MATH = Var.create_safe(
    "remarkMath", _var_is_local=False, _var_is_string=False)
_REMARK_GFM = Var.create_safe(
    "remarkGfm", _var_is_local=False, _var_is_string=False)
_REMARK_UNWRAP_IMAGES = Var.create_safe(
    "remarkUnwrapImages", _var_is_local=False, _var_is_string=False
)
_REMARK_TOC = Var.create_safe(
    "remarkToc",  _var_is_local=False, _var_is_string=False)

# Special rehype plugins.
_REHYPE_KATEX = Var.create_safe(
    "rehypeKatex", _var_is_local=False, _var_is_string=False
)
_REHYPE_RAW = Var.create_safe(
    "rehypeRaw", _var_is_local=False, _var_is_string=False)

_REHYPE_SLUG = Var.create_safe(
    "rehypeSlug", _var_is_local=False, _var_is_string=False)

_REHYPE_AUTOLINK_HEADINGS = Var.create_safe(
    "rehypeAutolinkHeadings", _var_is_local=False, _var_is_string=False)

_MOCK_ARG = Var.create_safe("", _var_is_string=False)

_REMARK_PLUGINS = Var.create_safe(
    [_REMARK_MATH, _REMARK_GFM, _REMARK_UNWRAP_IMAGES, _REMARK_TOC])

_REHYPE_PLUGINS = Var.create_safe(
    [_REHYPE_KATEX, _REHYPE_RAW, _REHYPE_SLUG, _REHYPE_AUTOLINK_HEADINGS])


class MyMarkdown(Markdown):
    def add_imports(self) -> ImportDict | list[ImportDict]:
        """Add imports for the markdown component.

        Returns:
            The imports for the markdown component.
        """
        from reflex.components.datadisplay.code import CodeBlock
        from reflex.components.radix.themes.typography.code import Code

        return [
            {
                "": "katex/dist/katex.min.css",
                "remark-math@5.1.1": ImportVar(
                    tag=_REMARK_MATH._var_name, is_default=True
                ),
                "remark-gfm@3.0.1": ImportVar(
                    tag=_REMARK_GFM._var_name, is_default=True
                ),
                "remark-unwrap-images@4.0.0": ImportVar(
                    tag=_REMARK_UNWRAP_IMAGES._var_name, is_default=True
                ),
                "rehype-katex@6.0.3": ImportVar(
                    tag=_REHYPE_KATEX._var_name, is_default=True
                ),
                "rehype-raw@6.1.1": ImportVar(
                    tag=_REHYPE_RAW._var_name, is_default=True
                ),
                "remark-toc@9.0.0": ImportVar(
                    tag=_REMARK_TOC._var_name, is_default=True
                ),
                "rehype-slug@6.0.0": ImportVar(
                    tag=_REHYPE_SLUG._var_name, is_default=True
                ),
                "rehype-autolink-headings@7.1.0": ImportVar(
                    tag=_REHYPE_AUTOLINK_HEADINGS._var_name, is_default=True
                ),
            },
            *[
                component(_MOCK_ARG)._get_imports()  # type: ignore
                for component in self.component_map.values()
            ],
            CodeBlock.create(theme="light")._get_imports(),  # type: ignore,
            Code.create()._get_imports(),  # type: ignore,
        ]

    def _render(self) -> Tag:
        tag = (
            Component
            ._render(self)
            .add_props(
                remark_plugins=_REMARK_PLUGINS,
                rehype_plugins=_REHYPE_PLUGINS,
            )
            .remove_props("componentMap", "componentMapHash")
        )
        tag.special_props.add(
            Var.create_safe(
                f"components={{{self._get_component_map_name()}()}}",
                _var_is_local=True,
                _var_is_string=False,
            ),
        )
        return tag
```

这样，我们按照官方的写法重新写一下这两个函数即可实现插件的引用了！还是比较方便的。

## 如何自动生成TOC

只需要在你的Markdown文件的特定位置写上`## Contents`，它就会自动帮你生成了，非常方便。
