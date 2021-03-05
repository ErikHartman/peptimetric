
import os
import sys

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output

app = dash.Dash(external_stylesheets=[dbc.themes.SANDSTONE])

NAVBAR_BRAND_STYLE = {
    "font-size": "20px",
    "margin-left": "5px" "!important"
 }
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": "4rem",
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "light",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}
navbar = dbc.NavbarSimple(
    children=[
        
        dbc.NavItem(dbc.NavLink("Add files", href="#")),
        dbc.NavItem(dbc.NavLink("Documentation", href="#")),
        dbc.NavItem(dbc.NavLink("FAQ", href="#"))
    ],
    brand="PeptiDiff",
    brand_style=NAVBAR_BRAND_STYLE,
    brand_href="#",
    color="light",
    dark=False,
)

sidebar = html.Div(
    [
        html.H5("Visualization tool for proteomics and peptidomics", className="display-8"),
        html.Hr(),
        dbc.Nav([
            dbc.NavLink("Menu 1", href="/", active="exact"),
            dbc.NavLink("Menu 2", href="/menu-2", active="exact"),
            dbc.NavLink("Menu 3", href="/menu-3", active="exact"),
        ],
        vertical=True,
        pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)
"""
#TODO: Fix card layout

card = dbc.Card(
    dbc.CardBody(
        [
            html.H4("Help", className="card-title"),
            html.H6("Using PeptiDiff", className="card-subtitle"),
            html.P(
                "Load files ",
                className="card-text",
            ),
        ]
    ),
    style={"width": "18rem", "height":"30rem",  "align":'center'}, color='light'
)
"""
content = html.Div(id="page-content", children=[], style=CONTENT_STYLE)

app.layout = html.Div([
    dcc.Location(id="url"),
    navbar,
    sidebar,
    content,
    
])
@app.callback(
    Output("page-content", "children"),
    [Input("url", "pathname")]
)
def render_menu_content(pathname):
    if pathname == "/":
        return [
                html.H1('Menu 1',
                        style={'textAlign':'center'}),
                ]
    elif pathname == "/menu-2":
        return [
                html.H1('Menu 2',
                        style={'textAlign':'center'}),
                ]
    elif pathname == "/menu-3":
        return [
                html.H1('Menu 3',
                        style={'textAlign':'center'}),
        ]

if __name__ == '__main__':
    app.run_server(debug=True)
