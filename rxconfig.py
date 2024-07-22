import reflex as rx

config = rx.Config(
    app_name="SequenceTool",
    tailwind={
        "theme": {
            "extend": {},
        },
        "plugins": ["@tailwindcss/typography"],
    },
)
