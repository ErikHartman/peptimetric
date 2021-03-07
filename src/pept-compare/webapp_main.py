
import os
import sys

import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output

app = dash.Dash(external_stylesheets=[dbc.themes.SANDSTONE])

NAVBAR_STYLE = {
    "left": 0,
    "background-color": "light",
    
    }

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": "2rem",
    "left": 0,
    "bottom": 0,
    "width": "12rem",
    "padding": "2rem 1rem",
    "background-color": "light",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}
navbar = dbc.Navbar(
    [
        dbc.NavbarBrand("Eriks och Simons kandidatarbete"),
            dbc.Nav(
            [
                dbc.DropdownMenu(
                    label="Files",
                    children=[
                        dbc.DropdownMenuItem("View files", href="/View-files"),
                        dbc.DropdownMenuItem("Add files", href="/Add-files"),
                    ],
                    className="mr-1",
                    nav=True,
                    in_navbar=True,
                ),
                dbc.Button(
                    "Export data",
                    className="mr-1",
                    color='info',
                    href="/Export-data"
                ),
                dbc.NavLink("Documentation", href="/Documentation", active="exact"),
                dbc.NavLink("FAQ", href="/FAQ", active="exact"),
                dbc.NavLink("Feedback", href="/Feedback", active="exact"),
            ],
            navbar=True,
            className="ml-auto",
            style=NAVBAR_STYLE
        ),   
    ]
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
        className="mb-4"
        ),
    ],
    style=SIDEBAR_STYLE,
)

search_input = html.Div([
    html.Div([
        dbc.Input(
            id='search_input',
            placeholder='Search protein...',
            debounce=True,
            minLength=0, maxLength=30,
            size = '20',
        )

    ])
])

how_to_use_card = dbc.Card(
    dbc.CardBody(
        [
            html.H4("How to use", className="card-title mb-4 font-weight-bold"),
            html.Hr(),
            html.H6("Using PeptiDiff", className="card-subtitle"),
            html.P(
                "Load files ",
                className="card-text",
            ),
        ]
    ),
    style={'height':'22rem'},
    color='white', className="border-dark mb-4"
)

protein_fig = dbc.Card(
    dbc.CardBody(
        [
            html.H4("Protein View", className="card-title mb-4 font-weight-bold"),
            html.Hr(),
        ]
    ),
    style={'height':'22rem'},
    color='white', className="border-dark"
)


peptide_fig = dbc.Card(
    dbc.CardBody(
        [
            html.H4("Peptide View", className="card-title mb-4 font-weight-bold"),
            html.Hr(),
        ]
    ),
    style={'height':'17rem'},
    color='white', className="border-dark"
)

columns = ['Peptide','Start','End','Intensity 1',' Intensity 2']
protein_info = dash_table.DataTable(
        id='protein_info',
        columns=[{"name": str(i), "id": str(i)} for i in columns],
        style_header={
        'fontWeight': 'bold',
        'backgroundColor':'primary',
    },
    )

peptide_info = dash_table.DataTable(
        id='peptide_info',
        columns=[{"name": str(i), "id": str(i)} for i in columns],
        style_header={
        'fontWeight': 'bold'
    },
    ),

content = html.Div(id="page-content", children=[], style=CONTENT_STYLE)

app.layout = dbc.Container([
    dcc.Location(id="url"),
    #Navbar
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4"),
    ]),
    #Search input 
    dbc.Row([
        dbc.Col(search_input, width={"size":0, "offset":8})
    ]),
    #Protein and how to row
    dbc.Row([
        dbc.Col(sidebar, width={'size':0}),
        dbc.Col(protein_fig, width={'size':5}),
        dbc.Col(protein_info, width={'size':1}),
        dbc.Col(how_to_use_card , width={'size':3, "offset":3}),
    ]),
    #Peptide row
    dbc.Row([
        dbc.Col(peptide_fig, width={'size': 8}),
        dbc.Col(peptide_info, width={'size':2})
    ])
    
    
])


if __name__ == '__main__':
    app.run_server(debug=True)
