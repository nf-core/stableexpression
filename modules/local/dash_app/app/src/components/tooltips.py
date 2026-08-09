import dash_mantine_components as dmc


def get_tooltip(
    classname: str, label: str, position: str = "bottom", multiline: bool = True
):
    return dmc.Tooltip(
        target=f".{classname}",
        label=label,
        multiline=multiline,
        position=position,
        color="grey",
        withArrow=True,
        arrowSize=8,
        zIndex=20000,
        radius=4,
        transitionProps={
            "transition": "fade",
            "duration": 200,
            "timingFunction": "ease",
        },
    )


genes_tabitem_tooltip = get_tooltip(
    classname="genes-tabitem", label="Distribution of normalised counts gene per gene"
)

samples_tabitem_tooltip = get_tooltip(
    classname="samples-tabitem",
    label="Distribution of normalised counts sample per sample",
)

settings_button_tooltip = get_tooltip(
    classname="settings-button",
    label="Open settings to select genes / samples and to customise display",
)

tooltips_to_load = [
    genes_tabitem_tooltip,
    samples_tabitem_tooltip,
    settings_button_tooltip,
]
