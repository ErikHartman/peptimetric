import os
import sys
import base64
import datetime
import io
import json
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State

from methods import make_peptide_dfs, concatenate_dataframes, merge_dataframes, apply_cut_off, create_protein_list, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly
from methods import amino_acid_piecharts, common_family
app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(id='page-content', fluid=True, className='vh-100'),
 ])

#---------------------------------------PAGE-ELEMENTS------------------------------------------------

modal_file = html.Div([
    dbc.Button("Files", id="open-modal-file", color='info', className="mr-1"),
        dbc.Modal([
                dbc.ModalHeader("Files", className="font-weight-bold"),
                    dbc.Row([
                        dbc.Col(dbc.ModalBody('Group 1', className='ml-auto text-center')),
                        dbc.Col(dbc.ModalBody('Group 2', className='ml-auto text-center')),
                    ]),
                    dbc.Row([
                        dbc.Col(dcc.Upload(id='upload-data-1', children=dbc.Button('Select files'), multiple=True), className="text-center"),
                        dbc.Col(dcc.Upload(id='upload-data-2', children=dbc.Button('Select files'), multiple=True, ), className="text-center"),
                        
                    ]),    
                    dbc.Row([
                        dbc.Col(dbc.ModalBody(id='output-filename-1', className='ml-auto text-center')),
                        dbc.Col(dbc.ModalBody(id='output-filename-2', className='ml-auto text-center')),
                    ]),
                dbc.ModalFooter(
                    dbc.Button("Generate protein graph", color = 'primary', id="close-modal-file", className="ml-auto", n_clicks_timestamp=0)
                ),
            ],
            id="modal-file",
            centered=True,
            scrollable=True,
        )])

modal_export_data = html.Div([
    dbc.Button("Export Data", id="open-modal-export-data", color='info', className="mr-1"),
    dbc.Modal([
                dbc.ModalHeader("Export data", className="font-weight-bold"),

                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-export-data", className="ml-auto")
                ),
            ],
            id="modal-export-data",
            centered=True,
              )
])

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
        modal_export_data,
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
        ),   
    ],
    
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
)

amino_acid_pie_dropdown = dbc.DropdownMenu(label='View', color="light",
children= [
    dbc.DropdownMenuItem('Complete proteome', id='complete-proteome', n_clicks_timestamp=0),
    dbc.DropdownMenuItem('Selected protein', id='selected-protein', n_clicks_timestamp=0)
]
)

search_protein = html.Div([
    html.Div([
        dbc.Input(
            id='search-protein',
            placeholder='Search protein...',
            name = 'text',
            debounce=True,
            inputMode='latin',
            minLength=0, maxLength=30,
            size = '20',
            list = 'protein-list',
            className="ml-auto",
        )
    ]),
    html.Datalist( id = 'protein-list', children=[])
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

protein_fig = html.Div([
        html.H3('Protein View'),
        dbc.Row([
            dbc.Col(search_protein),
            dbc.Col(dbc.Button('Show protein families', id='show-pfam', className='ml-auto', n_clicks_timestamp=0)),
            dbc.Col(dbc.Button('Show error bars', id='show-protein-error-bars', className='ml-auo', n_clicks_timestamp=0)),
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='protein-fig', figure={})
        )
        
        ])


peptide_fig = html.Div([
        html.H3('Peptide View'),
        dbc.Row([
            dbc.Col(dbc.Button('Show weight', id='show-peptide-weight', className='ml-auto')),
            dbc.Col(dbc.Button('Show difference trace', id='show-difference-trace', className='ml-auo')),
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='peptide-fig', figure={})
        )
        ])

amino_acid_figs = html.Div([
        html.H3('Amino Acid Profile'),
        amino_acid_pie_dropdown,
        html.P(id='complete-or-selected'),
        dcc.Loading(type='cube', color = '#76b382',
            children=[ dbc.Row([
                dbc.Col([
                    dcc.Graph(id='complete-aa-seq-fig-g1', figure={})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='first-aa-fig-g1', figure={})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='last-aa-fig-g1', figure={})], width={'size':3}),
            ]),
            dbc.Row([
                dbc.Col([
                    dcc.Graph(id='complete-aa-seq-fig-g2', figure={})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='first-aa-fig-g2', figure={})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='last-aa-fig-g2', figure={})], width={'size':3}),
            ])]
        )
        ])
    
table_header = [html.Thead(html.Tr([html.Th("Protein info")]))]
row1 = html.Tr([html.Td("Number of proteins"), html.Td("0", id='number-of-proteins')])
row2 = html.Tr([html.Td("Average number of peptides per protein"), html.Td("0", id='average-nbr-peptides-per-protein')])
row3 = html.Tr([html.Td("Average peptide length"), html.Td("0", id='average-peptide-length')])
table_body = [html.Tbody([row1, row2, row3])]

protein_info = dbc.Table(table_header + table_body, 
    bordered=True,
    hover=True,
    responsive=True,
    striped=True,)    



peptide_info = html.Div(
   id='peptide-info-table'),

#---------------------------PAGES---------------------------------------------------------------
main_page = dbc.Container([
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4")
    ]),
    dbc.Row([
        dbc.Col(how_to_use_collapse , width={'size':12}),
    ]),
    dbc.Row([
        dbc.Col(protein_fig, width={'size':8}),
        dbc.Col(protein_info, width={'size':4}),
    ]),
    dbc.Row([
        dbc.Col(peptide_fig, width={'size': 8}),
        dbc.Col(peptide_info, width={'size':4})
    ]),

    dbc.Row([
        dbc.Col(amino_acid_figs),
    ])      
], fluid=True)

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

df_g = []
def update_file_list(contents, filename):
    s = ""
    i=1
    if filename:
        for f in filename:
            s = s + (f"S{i}:{f}" + "\n")
            i+=1  
        master_df = update_data_frame(contents, filename)
        df_g.append(master_df)
    return s

def update_data_frame(contents, filename):
    decoded_list = []
    for f, content in zip(filename, contents):
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        decoded_list.append(io.BytesIO(decoded))
    dfs = make_peptide_dfs(decoded_list)
    master_df = concatenate_dataframes(dfs)
    return master_df

protein_lists = []
def create_protein_fig(n_clicks, n_clicks_pfam, n_clicks_stdev):
    triv_names = []
    protein_fig = {}
    protein_list = []
    protein_df = pd.DataFrame()
    if n_clicks:
        g1 = df_g[-2]
        g2 = df_g[-1]
        master = merge_dataframes(g1,g2)
        protein_list = create_protein_list(master)
        protein_lists.append(protein_list)
        protein_list_cutoff = apply_cut_off(protein_list, nbr_of_peptides=5, area=1000000, spectral_count=4)
        protein_fig, protein_df = protein_graphic_plotly(protein_list_cutoff, difference_metric='area_sum')
        if len(protein_list_cutoff) > 1:
            for protein in protein_list_cutoff:
                triv_names.append(html.Option(value=protein.get_trivial_name()))
    if n_clicks > n_clicks_pfam and n_clicks > n_clicks_stdev and df_g and len(df_g) >= 2:
        return protein_fig, triv_names
    elif n_clicks_pfam or n_clicks_stdev and len(df_g) >= 2:
        return show_pfam_or_stdev(n_clicks_pfam, n_clicks_stdev, protein_list, protein_fig, protein_df) , triv_names
    else:
        return {}, []

peptide_lists=[]
def create_peptide_fig(clickData, search_protein):
    protein_accession = ''
    search_text = ''
    if search_protein != '' and len(protein_lists) > 0:
        for protein in protein_lists[-1]:
            if search_protein == protein.get_trivial_name():
                protein_accession = protein.get_id()
                search_text = ''
    elif clickData:
        protein_accession = clickData['points'][0]['customdata'][-1]
        search_text=''

    if protein_accession != '':
        peptide_list = create_peptide_list(protein_lists[-1], str(protein_accession))
        peptide_lists.append(peptide_list)
        peptide_fig = peptide_graphic_plotly(peptide_list)
        peptide_info_columns = ['Peptide','Start','End','Intensity_g1','Intensity_g2']
        df_peptide_info = pd.DataFrame(columns=peptide_info_columns)
        for peptide in peptide_list:
            df_peptide_info = df_peptide_info.append({'Peptide': str(peptide.get_sequence()), 'Start': peptide.get_start(),'End': peptide.get_end(), 'Intensity_g1': "{:.2e}".format(peptide.get_area()[0]), 
            'Intensity_g2': "{:.2e}".format(peptide.get_area()[1])}, ignore_index=True)
        print(df_peptide_info)
        df_peptide_info.sort_values(by=['Intensity_g1', 'Intensity_g2'], ascending=False, inplace=True)
        datatable = dash_table.DataTable(
            data = df_peptide_info.to_dict('rows'),
            columns=[{"name": str(i), "id": str(i)} for i in df_peptide_info.columns],
            sort_action='native',
            fixed_rows={'headers': True},
            filter_action='native',
            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : 'rgb(182, 224, 194)'
            }
            ],
            style_header={
                'textAlign':'center',
                'fontWeight': 'bold',
            },
            style_cell={
                'textAlign':'left',
                'padding':'5px',
                'maxWidth': 95,
            },
            style_table={'height': '400px', 'width':'100%', 'overflowY': 'auto'}
    )   
        return peptide_fig, search_text, html.Div([datatable])
    else:
        return {}, search_text, []


def amino_acid_dropdown(n_clicks_complete_proteome, n_clicks_selected_protein):
    if n_clicks_complete_proteome > n_clicks_selected_protein and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(protein_lists[-1], peptide_or_protein_list = 'protein_list')
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Complete proteome'
    if n_clicks_selected_protein > n_clicks_complete_proteome and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_lists[-1], peptide_or_protein_list = 'peptide_list')
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Selected Protein'
    else:
        return  {},{},{}, {},{},{}, ''


def show_pfam_or_stdev(n_clicks_pfam, n_clicks_stdev, protein_list, fig, df):
    if n_clicks_pfam > n_clicks_stdev:
        for p1 in protein_list:
            for p2 in protein_list:
                if common_family(p1.get_protein_family(), p2.get_protein_family())[0]:
                    x0 = p1.get_area_sum()[2]
                    x1 = p2.get_area_sum()[2]
                    y0 = p1.get_area_sum()[0]
                    y1 = p2.get_area_sum()[0]
                    fig.add_shape(type="line",x0=x0, y0=y0, x1=x1, y1=y1, line=dict(color="firebrick",width=1, dash='dash'))
        return fig
    elif n_clicks_stdev and len(protein_list) > 2:
        fig.update_traces(error_x= dict(array=df['g2_stdev'].array, thickness=1), error_y=dict(array=df['g1_stdev'].array, thickness=1))
        return fig
    else:
        return {}



app.callback(
    Output('complete-aa-seq-fig-g1', 'figure'),
    Output('first-aa-fig-g1', 'figure'),
    Output('last-aa-fig-g1', 'figure'),
    Output('complete-aa-seq-fig-g2', 'figure'),
    Output('first-aa-fig-g2', 'figure'),
    Output('last-aa-fig-g2', 'figure'),
    Output('complete-or-selected', 'children'),
    Input('complete-proteome', 'n_clicks_timestamp'),
    Input('selected-protein','n_clicks_timestamp')
)(amino_acid_dropdown)

app.callback(
    Output('peptide-fig', 'figure'),
    Output('search-protein', 'value'),
    Output('peptide-info-table', 'children'),
    Input('protein-fig', 'clickData'),
    Input('search-protein', 'value'),
)(create_peptide_fig)

app.callback(
    Output("output-filename-1", "children"),
    [Input("upload-data-1", "contents"), Input("upload-data-1", "filename")],
)(update_file_list)

app.callback(
    Output("output-filename-2", "children"),
    [Input("upload-data-2", "contents"), Input("upload-data-2", "filename")],
)(update_file_list)

app.callback(
    Output('protein-fig', 'figure'),
    Output('protein-list', 'children'),
    [Input("close-modal-file", "n_clicks_timestamp"),
    Input('show-pfam', 'n_clicks_timestamp'),
    Input('show-protein-error-bars','n_clicks_timestamp')],
    )(create_protein_fig)

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
    Output("modal-export-data", "is_open"),
    [Input("open-modal-export-data", "n_clicks"), Input("close-modal-export-data", "n_clicks")],
    [State("modal-export-data", "is_open")],
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
