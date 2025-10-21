import dash_mantine_components as dmc

from src.utils import style

from src.components.settings import genes, samples

drawer = dmc.Drawer(
    children=[
        genes.sidebar_stack,
        samples.sidebar_stack,
    ],
    id="drawer",
    opened=False,
    position="right",
    withCloseButton=True,
    closeOnEscape=True,
    overlayProps=dict(backgroundOpacity=0),
    trapFocus=False,
    zIndex=10000,
    style=style.SIDEBAR,
)
