from dash_extensions.enrich import Input, Trigger, Output, State, callback


##############################################
##############################################
# CALLBACKS
##############################################
##############################################


def register_callbacks():
    @callback(
        Output("drawer", "opened"),
        Trigger("settings-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def open_drawer():
        return True

    @callback(
        Output("sidebar-genes-items", "style"),
        Output("sidebar-samples-items", "style"),
        Input("tabs", "value"),
        State("sidebar-genes-items", "style"),
        State("sidebar-samples-items", "style"),
    )
    def manage_drawer_content(
        tabs_value: str, gene_stack_style: dict, sample_stack_style: dict
    ):
        if tabs_value == "genes":
            gene_stack_style["display"] = "block"
            sample_stack_style["display"] = "none"
        else:  # tabs_value ==  'samples':
            gene_stack_style["display"] = "none"
            sample_stack_style["display"] = "block"
        return gene_stack_style, sample_stack_style
