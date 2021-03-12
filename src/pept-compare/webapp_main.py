
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

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)


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
                dbc.ModalHeader("Files", className="font-weight-bold"),
                dbc.Row([
                    dbc.Col(dbc.ModalBody('Group 1')),
                    dbc.Col(dbc.ModalBody('Group 2')),
                ]),
                
                dbc.Row([
                    dbc.Col([
                    dcc.Upload(
                    id='upload-data',
                    children=dbc.Button('Select files', className="ml-auto"),
                    multiple=True)
                    ]),
                    dbc.Col([
                    dcc.Upload(
                    id='upload-data',
                    children=dbc.Button('Select files', className="ml-auto"),
                    multiple=True),
                    ]),
                dbc.Row([
                    dbc.Col([
                        html.Div(id='output-data-upload')
                    ])
                ])
                ]),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-file", className="ml-auto")
                ),
            ],
            id="modal-file",
            centered=True,
        )])


modal_other_graphics = dbc.Modal([
                dbc.ModalHeader("Other Graphics", className="font-weight-bold"),

                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-other-graphics", className="ml-auto")
                ),
            ],
            id="modal-other-graphics",
            centered=True,
              )

modal_view_setings = dbc.Modal([
                dbc.ModalHeader("View Settings", className="font-weight-bold"),

                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-view-settings", className="ml-auto")
                ),
            ],
            id="modal-view-settings",
            centered=True,
              )

modal_cutoff = dbc.Modal([
                dbc.ModalHeader("Cutoff settings", className="font-weight-bold"),

                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-cutoff", className="ml-auto")
                ),
            ],
            id="modal-cutoff",
            centered=True,
              )

navbar = dbc.Navbar(
    [
        dbc.NavbarBrand("Eriks och Simons kandidatarbete"),
        modal_file,
        dbc.Button("Export data", className="mr-1", color='info', href="/Export-data"),
        dbc.DropdownMenu(label="Settings",
            children=[
                dbc.DropdownMenuItem("View", id="open-modal-view-settings"),
                modal_view_setings,
                dbc.DropdownMenuItem("Cutoffs", id="open-modal-cutoff"),
                modal_cutoff,
                dbc.DropdownMenuItem("Other graphics", id="open-modal-other-graphics"),
                modal_other_graphics,
            ]
        ),
            dbc.Nav(
            [
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
how_to_use_collapse = html.Div(
    [
        dbc.Button(
            "How to Use",
            id="collapse-button",
            className="mb-3",
            color="info",
        ),
        dbc.Collapse(
            dbc.Card([
            dbc.CardHeader("How to use!"),
            dbc.CardBody("This will display steps and links on how to use the app."),
            ]),
            id="collapse",
            className="mb-4"
        ),
    ]
)

protein_fig = dbc.Card(
    dbc.CardBody(
        [
            html.H4("Protein View", className="card-title mb-4 font-weight-bold"),
            html.Hr(),
        ]
    ),
    style={'height':'34rem'},
    color='white', className="border-dark mb-4"
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
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4"),
    ]),
    dbc.Row([
        dbc.Col(how_to_use_collapse , width={'size':12, "offset":0}),
    ]),
    dbc.Row([
        dbc.Col(search_input, width={"size":2, "offset":3}, className="mb-4"),
    ]),
    dbc.Row([
        dbc.Col(protein_fig, width={'size':8}),
        dbc.Col(protein_info, width={'size':1}),
    ]),
    dbc.Row([
        dbc.Col(peptide_fig, width={'size': 8}),
        dbc.Col(peptide_info, width={'size':2})
    ])   
])

FAQ_page = html.Div([
    navbar,
    html.H1('FAQ'),
    html.Br(),
    dcc.Link('Go back to home', href='/')
])

documentation_page = html.Div([
    navbar,
    html.H1('Documentation'),
    html.Br(),
    dcc.Link('Go back to home', href='/')
])

#-----------------DEFS AND CALLBACKS--------------------------------------------------------------

def display_page(pathname):
    if pathname == '/FAQ':
        return FAQ_page
    elif pathname == '/Documentation':
        return documentation_page
    else: 
        return main_page


def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

def parse_data(contents, filename):
    try:
        if 'csv' in filename:
            df = pd.read_csv(filename)
        elif 'xls' in filename:
            df = pd.read_excel(filename)
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    print(filename)
    return html.Div([
        html.H5(filename)])
    
@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n) for c, n in
            zip(list_of_contents, list_of_names)]
        return children

app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])(display_page)

app.callback(
    Output("modal-file", "is_open"),
    [Input("open-modal-file", "n_clicks"), Input("close-modal-file", "n_clicks")],
    [State("modal-file", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-other-graphics", "is_open"),
    [Input("open-modal-other-graphics", "n_clicks"), Input("close-modal-other-graphics", "n_clicks")],
    [State("modal-other-graphics", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-cutoff", "is_open"),
    [Input("open-modal-cutoff", "n_clicks"), Input("close-modal-cutoff", "n_clicks")],
    [State("modal-cutoff", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-view-settings", "is_open"),
    [Input("open-modal-view-settings", "n_clicks"), Input("close-modal-view-settings", "n_clicks")],
    [State("modal-view-settings", "is_open")],
)(toggle_modal)

app.callback(
    Output("collapse", "is_open"),
    [Input("collapse-button", "n_clicks")],
    [State("collapse", "is_open")],
)(toggle_collapse)



if __name__ == '__main__':
    app.run_server(debug=True)
