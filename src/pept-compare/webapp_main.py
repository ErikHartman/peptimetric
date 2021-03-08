
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
from dash.dependencies import Input, Output, State

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE])


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
 ])


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
#---------------------------------------PAGE-ELEMENTS------------------------------------------------
modal_file = html.Div([
        dbc.Button("Files", id="open-modal-file", color='info', className="mr-1"),
        dbc.Modal([
                dbc.ModalHeader("Files"),
                dbc.Row([
                    dbc.Col(dbc.ModalBody('Group 1')),
                    dbc.Col(dbc.ModalBody('Group 2')),
                ]),
                dbc.Row([
                dbc.Col(
                    dcc.Upload(id="upload-group-1",
                    children=html.Button('Upload Files'),
                    className="ml-auto",
                    ),
                ),
                dbc.Col(
                    dcc.Upload(id="upload-group-2",
                    children=html.Button('Upload Files'),
                    className="ml-auto",),
                    ),
                ]),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-file", className="ml-auto")
                ),
            ],
            id="modal-file",
            is_open=False,
        )
])
navbar = dbc.Navbar(
    [
        dbc.NavbarBrand("Eriks och Simons kandidatarbete"),
            dbc.Nav(
            [
                modal_file,
                dbc.Button("Export data", className="mr-1", color='info', href="/Export-data"),
                dbc.NavLink("Home", href="/", active="exact"),
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

#---------------------------PAGES---------------------------------------------------------------
main_page = dbc.Container([
    dcc.Location(id="url"),
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4"),
    ]),
    dbc.Row([
        dbc.Col(search_input, width={"size":0, "offset":2}, className="mb-4"),
    ]),
    dbc.Row([
        dbc.Col(sidebar, width={'size':0}),
        dbc.Col(protein_fig, width={'size':5}),
        dbc.Col(protein_info, width={'size':1}),
        dbc.Col(how_to_use_card , width={'size':3, "offset":3}),
    ]),
    dbc.Row([
        dbc.Col(peptide_fig, width={'size': 8}),
        dbc.Col(peptide_info, width={'size':2})
    ])   
])

FAQ_page = html.Div([
    dcc.Location(id="url"),
    navbar,
    html.H1('FAQ'),
    html.Br(),
    dcc.Link('Go back to home', href='/')
])

documentation_page = html.Div([
    dcc.Location(id="url"),
    navbar,
    html.H1('Documentation'),
    html.Br(),
    dcc.Link('Go back to home', href='/')
])

@app.callback(dash.dependencies.Output('page-content', 'children'),
              [dash.dependencies.Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/FAQ':
        return FAQ_page
    elif pathname == '/Documentation':
        return documentation_page
    elif pathname == '/':
        return main_page
    else: 
        return main_page

@app.callback(
    Output("modal-file", "is_open"),
    [Input("open-modal-file", "n_clicks"), Input("close-modal-file", "n_clicks")],
    [State("modal-file", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


if __name__ == '__main__':
    app.run_server(debug=True)
