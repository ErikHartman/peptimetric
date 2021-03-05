
import pandas as pd
import plotly.express as px 
import plotly.graph_objects as go

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import sys, os

app = dash.Dash(external_stylesheets=[dbc.themes.SANDSTONE])

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f8f8",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

sidebar = html.Div(
    [
        html.H2("PeptiDiff", className="display-4"),
        html.Hr(),
        html.H4("Visualization tool for proteomics and peptidomics", className="display-8"),
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

content = html.Div(id="page-content", children=[], style=CONTENT_STYLE)

app.layout = html.Div([
    dcc.Location(id="url"),
    sidebar,
    content,
    

    dcc.Textarea(
        id='textarea-example',
        placeholder='Serach protein...', 
        value='Hello World',
        style={'width': '40%', 'height': 500, 'font_family': 'Roboto'} , 
   ),
    html.Div(id='textarea-example-output', style={'whiteSpace': 'pre-line'})
    
])

if __name__ == '__main__':
    app.run_server(debug=True)
