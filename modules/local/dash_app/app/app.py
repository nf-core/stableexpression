import os
import dash_mantine_components as dmc
from dotenv import load_dotenv

from dash_extensions.enrich import (
    DashProxy,
    html,
    ServersideOutputTransform,
    TriggerTransform,
)
from dash_extensions.logging import NotificationsLogHandler

from src.utils import config, style
from src.components import stores, tooltips
from src.components import top, right_sidebar
from src.callbacks import common, genes, samples

load_dotenv("./.env")
debug = True if os.getenv("DEBUG") is not None else False


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
    suppress_callback_exceptions=(not debug),
    update_title=config.UPDATE_TITLE,
    external_stylesheets=[dmc.styles.ALL],
    transforms=[TriggerTransform(), ServersideOutputTransform()],
)

# -------------------- LAYOUT --------------------


def serve_layout():
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
            )
        ]
    )


app.layout = serve_layout

# -------------------- IMPORTING CALLBACKS --------------------

common.register_callbacks()
genes.register_callbacks()
samples.register_callbacks()


if __name__ == "__main__":
    logger.info("Running server")
    # setting prune_errors to False avoids error message pruning
    # in order to get original tracebacks
    # (very useful for debugging)
    prune_errors = False if debug else True
    app.run(
        debug=debug,
        host=config.HOST,
        port=config.PLOTLY_APP_PORT,
        dev_tools_prune_errors=prune_errors,
    )
