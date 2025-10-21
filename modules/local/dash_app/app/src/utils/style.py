LAYOUT = {
    "left": "0px",
    "top": "0px",
    "position": "absolute",
    "width": "100%",
    "height": "100%",
}

HEADER_HEIGHT = "5em"

TAB = {
    "position": "fixed",
    "top": 0,
    "left": 10,
    "right": 0,
    "width": "100%",
    "height": "100%",
    #'zIndex': '1001',
}

HEADER_TABLIST = {
    "position": "fixed",
    "top": 0,
    "left": 10,
    "right": 0,
    "width": "60%",
    "height": HEADER_HEIGHT,
    #'zIndex': '1001'
}

HEADER_TABLIST_ITEM = {
    #'width': '15vh',
    # "text-align": "center",
    "paddingRight": "20px",
    #'paddingTop': '26px',
    #'paddingBottom': '26px',
    #'width': LEFT_SIDEBAR_WIDTH
}

TABS_PANEL = {"margin-top": HEADER_HEIGHT, "height": f"calc(100% - {HEADER_HEIGHT})"}


SETTINGS_BUTTON = {
    "right": "20px",
}

SIDEBAR_WIDTH = "15em"

SIDEBAR = {
    "position": "fixed",
    "top": HEADER_HEIGHT,
    "bottom": 0,
    "width": SIDEBAR_WIDTH,
    "height": "100vh",
    "alignItems": "center",
}


DROPDOWN = {"marginTop": "10px", "paddingLeft": "4.2em", "paddingRight": "4.5em"}

STACK_SUBSECTION_TITLE = {"marginBottom": "-20px"}

AG_GRID = {
    "height": "calc(100% - 10px)",
    "top": HEADER_HEIGHT,
    "paddingTop": "10px",
    "marginRight": "15px",
    "paddingRight": "15px",
    "marginLeft": "5px",
}

GRAPH = {
    #'width': '100vh',
    "top": HEADER_HEIGHT,
    "marginLeft": "0px",
    "marginRight": "3em",
    "marginTop": "2px",
    "marginBottom": "3px",
    "display": "none",
}
