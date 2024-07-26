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


component_map = {
    "table": lambda children: rx.table.root(children, variant="surface"),
    "thead": lambda children: rx.table.header(children),
    "tbody": lambda children: rx.table.body(children),
    "tr": lambda children: rx.table.row(children),
    "th": lambda text: rx.table.column_header_cell(text, justify="center"),
    "td": lambda text: rx.table.cell(text, justify="center"),
    "codeblock": lambda text, **props: rx.code_block(text, **props, theme="light", show_line_numbers=False, wrap_long_lines=False),
}

mymarkdown = MyMarkdown.create
