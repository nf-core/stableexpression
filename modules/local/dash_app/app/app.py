import socket

import dash_mantine_components as dmc
from dash_extensions.enrich import (
    DashProxy,
    ServersideOutputTransform,
    TriggerTransform,
    html,
)
from dash_extensions.logging import NotificationsLogHandler
from src.callbacks import common, genes, samples
from src.components import right_sidebar, stores, tooltips, top
from src.utils import config, style

DEBUG = True

# -------------------- SETUP LOGGING --------------------

log_handler = NotificationsLogHandler()
logger = log_handler.setup_logger(__name__)

# -------------------- APP --------------------
# init the application
logger.info("Creating app")

app = DashProxy(
    __name__,
    title=config.APP_TITLE,
    prevent_initial_callbacks="initial_duplicate",
    suppress_callback_exceptions=(not DEBUG),
    update_title=config.UPDATE_TITLE,
    external_stylesheets=[dmc.styles.ALL],
    transforms=[TriggerTransform(), ServersideOutputTransform()],
)

# -------------------- LAYOUT --------------------


def serve_layout() -> dmc.MantineProvider:
    """Return the layout of the app."""
    return dmc.MantineProvider(
        children=[
            html.Div(
                [
                    top.header,
                    right_sidebar.drawer,
                    *stores.stores_to_load,
                    *tooltips.tooltips_to_load,
                ]
                + log_handler.embed(),
                id="layout",
                style=style.LAYOUT,
            ),
        ],
    )


app.layout = serve_layout

# -------------------- IMPORTING CALLBACKS --------------------

common.register_callbacks()
genes.register_callbacks()
samples.register_callbacks()

# -------------------- LAUNCH SERVER --------------------


def find_port(port: int) -> int:
    """Find a port not in use starting at given port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("localhost", port)) == 0:
            return find_port(port=port + 1)
        return port


if __name__ == "__main__":
    logger.info("Running server")
    # setting prune_errors to False avoids error message pruning
    # in order to get original tracebacks
    # (very useful for debugging)
    prune_errors = False if DEBUG else True
    app.run(
        debug=DEBUG,
        host=config.HOST,
        port=find_port(port=config.PLOTLY_APP_PORT),
        dev_tools_prune_errors=prune_errors,
    )
