import dash_mantine_components as dmc
from src.components.settings import genes, samples
from src.utils import style

drawer = dmc.Drawer(
    children=[
        genes.sidebar_stack,
        samples.sidebar_stack,
    ],
    id="drawer",
    opened=True,
    position="right",
    withCloseButton=True,
    closeOnEscape=True,
    overlayProps=dict(backgroundOpacity=0),
    trapFocus=False,
    zIndex=10000,
    style=style.SIDEBAR,
)
