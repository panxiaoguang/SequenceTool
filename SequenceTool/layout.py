from typing import Callable
import reflex as rx
from .components import navbar, footer


def template(Contents: Callable[[], rx.Component]) -> rx.Component:
    return rx.box(
        navbar(),
        Contents(),
        footer(),
    )
