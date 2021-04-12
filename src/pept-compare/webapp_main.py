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

from methods import make_peptide_dfs, concatenate_dataframes, merge_dataframes, apply_cut_off, create_protein_list, protein_graphic_plotly, create_peptide_list, stacked_samples_peptide
from methods import amino_acid_piecharts, common_family, all_sample_bar_chart

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(id='page-content', fluid=True, className='vh-100'),
 ])

#---------------------------------------PAGE-ELEMENTS------------------------------------------------
file_columns = ['Sample', 'File']

modal_file = html.Div([
    dbc.Button("Files", id="open-modal-file", color='info', className="mr-1"),
        dbc.Modal([
                dbc.ModalHeader("Files", className="font-weight-bold"),
                    dbc.Row([
                        dbc.Col(dbc.ModalBody('Group 1', className='ml-auto text-center')),
                        dbc.Col(dbc.ModalBody('Group 2', className='ml-auto text-center')),
                        
                    ]),
                    dbc.Row([
                        dbc.Col(dcc.Upload(id='upload-data-1', children=dbc.Button('Select files'), multiple=True), className="text-center ml-auto"),
                        dbc.Col(dcc.Upload(id='upload-data-2', children=dbc.Button('Select files'), multiple=True), className="text-center ml-auto"),
                        
                    ]),    
                    dbc.Row([
                        dbc.Col(dash_table.DataTable(
                                id = 'output-filename-1',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                    'if' : {'row_index':'odd'},
                                    'backgroundColor' : 'rgb(182, 224, 194)'}
                                ],
                                style_header={
                                    'textAlign':'center',
                                    'fontWeight': 'bold',
                                },
                                style_cell={
                                    'textAlign':'center',
                                },
                            )),
                        dbc.Col(dash_table.DataTable(
                                id = 'output-filename-2',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                    'if' : {'row_index':'odd'},
                                    'backgroundColor' : 'rgb(182, 224, 194)'}
                                ],
                                style_header={
                                    'textAlign':'center',
                                    'fontWeight': 'bold',
                                },
                                style_cell={
                                    'textAlign':'center',
                                },
                            )),
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
protein_tab = dbc.Form([
                    dbc.FormGroup([
                    dbc.Label("Cutoff 1", className="mr-2"),
                    dbc.Input('Enter number...', placeholder='0', type='number', className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 2', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 3', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 4', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",)
                ],
                inline=True),

peptide_tab = dbc.Form([
                    dbc.FormGroup([
                    dbc.Label("Cutoff 5", className="mr-2"),
                    dbc.Input('Enter number...', placeholder='0', type='number', className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 6', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 7', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Cutoff 8', className='mr-2'),
                    dbc.Input('Enter number...', placeholder='0', type="number", className='ml-auto', min=0),
                    ], className="mr-3",)
                ],
                inline=True),

modal_cutoff = dbc.Modal([
                dbc.ModalHeader("Cutoff settings", className="font-weight-bold"),
                dbc.Tabs([
                    dbc.Tab(protein_tab, label='Protein'),
                    dbc.Tab(peptide_tab, label='Peptide')
                ]),
                dbc.ModalFooter(
                    dbc.Button("Apply", id="close-modal-cutoff", className="ml-auto")
                ),
            ],
            id="modal-cutoff",
            size='sm',
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
            id="how-to-use-collapse-button",
            className="mb-3",
            color="info",
        ),
        dbc.Collapse(
            dbc.Card([
            dbc.CardHeader("How to use!"),
            dbc.CardBody("This will display steps and links on how to use the app."),
            ]),
            id="how-to-use-collapse",
            className="mb-4"
        ),
    ]
)

sample_collapse = html.Div(
    [
        dbc.Button(
            "Samples",
            id="sample-collapse-button",
            className="mb-3",
            color="info",
        ),
        dbc.Collapse(
            dbc.Card([
            dbc.Row([
                dbc.Col('Group 1'),
                dbc.Col('Group 2')
                ]),
            dbc.Row([
            dbc.Col(dash_table.DataTable(
                                id = 'sample-collapse-1',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                    'if' : {'row_index':'odd'},
                                    'backgroundColor' : 'rgb(182, 224, 194)'}
                                ],
                                style_header={
                                    'textAlign':'center',
                                    'fontWeight': 'bold',
                                },
                                style_cell={
                                    'textAlign':'center',
                                },
                            )),
            dbc.Col(dash_table.DataTable(
                                id = 'sample-collapse-2',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                    'if' : {'row_index':'odd'},
                                    'backgroundColor' : 'rgb(182, 224, 194)'}
                                ],
                                style_header={
                                    'textAlign':'center',
                                    'fontWeight': 'bold',
                                },
                                style_cell={
                                    'textAlign':'center',
                                },
                            )),
                ]),
            ]),
            id="sample-collapse",
            className="mb-4"
        ),
    ]
)

protein_fig = html.Div([
        html.H3('Protein View'),
        dbc.Row([
            dbc.Col(search_protein, width={'size':3}),
            dbc.Col(dbc.Checklist(
                className='ml-auto',
                id='protein-checklist',
                switch=True,
                inline=True,
                options=[
                    {'label': 'Show standard deviation', 'value': 'show-stdev'},
                    {'label': 'Show protein families', 'value': 'show-pfam'}
                    ],
                )                 
        )]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='protein-fig', figure={})
        )
        
        ])

all_samples_protein_fig = html.Div([
    dcc.Graph(id='hover-all-protein-samples', figure={}, style={'height': 300, 'width':500}
)])

peptide_fig = html.Div([
        html.H3('Peptide View'),
        dbc.Row([
            dbc.Col(
            dbc.DropdownMenu(label='peptide dropdown',
            children = [
                dbc.DropdownMenuItem("Sum", id="peptide-dropdown-sum", n_clicks_timestamp=0),
                dbc.DropdownMenuItem("Mean", id="peptide-dropdown-mean", n_clicks_timestamp=0),
            ]))                 
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='peptide-fig', figure={})
        )
        ])

amino_acid_radioitems = html.Div([
        dbc.Label("Select difference metric"), 
        dbc.RadioItems(
            options=[
                {"label": "Area", "value": 'area'},
                {"label": "Spectral count", "value": 'spectral_count'}, 
            ], 
            value=1,
            id="aa-radioitems", 
            inline=True,

        )
])

amino_acid_figs = html.Div([
        html.H3('Amino Acid Profile'),
        dbc.Row([
            dbc.Col(amino_acid_pie_dropdown),
            dbc.Col(amino_acid_radioitems)

        ]),
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
    


protein_info = html.Div(id = 'protein-info-table'), 



peptide_info = html.Div([
    html.Div(id='peptide-info-table'),
])
#---------------------------PAGES---------------------------------------------------------------
main_page = dbc.Container([
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4")
    ]),
    dbc.Row([
        dbc.Col(how_to_use_collapse , width={'size':8}),
        dbc.Col(sample_collapse, width={'size':4}),
    ]),
    dbc.Row([
        dbc.Col(protein_fig, width={'size':8}),
        dbc.Col([
            dbc.Row(protein_info),
            dbc.Row(all_samples_protein_fig)], 
            width={'size':4}),
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
    if filename:
        file_list = []
        i=0
        for f in filename:
            s = [f"S{i}", f]
            file_list.append(s)
            i+=1  
        master_df = update_data_frame(contents, filename)
        df_g.append(master_df)
        df = pd.DataFrame(file_list, columns = ['Sample', 'File'])
        print(df)
        return df.to_dict('rows'), df.to_dict('rows')
    else:
        return [],[]

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
def create_protein_fig(n_clicks, checkbox_values):
    triv_names = []
    protein_fig = {}
    protein_list = []
    protein_list_cutoff = []
    if n_clicks and len(df_g) >= 2:
        g1 = df_g[-2]
        g2 = df_g[-1]
        master = merge_dataframes(g1,g2)
        protein_list = create_protein_list(master)
        protein_lists.append(protein_list)
        protein_list_cutoff = apply_cut_off(protein_list, nbr_of_peptides=5, area=1000000, spectral_count=4)
        protein_info_columns = ['Protein','UniProt id','#peptides g1','#peptides g2','Intensity_g1','Intensity_g2', 'Protein family']
        df_protein_info = pd.DataFrame(columns=protein_info_columns)
        for protein in protein_list:
            df_protein_info = df_protein_info.append({'Protein': str(protein.get_trivial_name()), 'UniProt id': protein.get_id(),'#peptides g1': protein.get_nbr_of_peptides()[0], '#peptides g2': protein.get_nbr_of_peptides()[1], 
            'Intensity_g1': "{:.2e}".format(protein.get_area_sum()[0]), 'Intensity_g2': "{:.2e}".format(protein.get_area_sum()[2]), 'Protein family':protein.get_protein_family()},  ignore_index=True)
        df_protein_info.sort_values(by=['Intensity_g1', 'Intensity_g2'], ascending=False, inplace=True)
        datatable = dash_table.DataTable(
            data = df_protein_info.to_dict('rows'),
            columns=[{"name": str(i), "id": str(i)} for i in df_protein_info.columns],
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
                'maxWidth': 105,
                'minWidth': 105,
            },
            style_table={'height': '200px', 'width':'500px', 'overflowY': 'auto','overflowX':'auto'}
    )   
        if len(protein_list_cutoff) > 1:
            for protein in protein_list_cutoff:
                triv_names.append(html.Option(value=protein.get_trivial_name()))
    if len(df_g) >= 2:
        if checkbox_values and 'show-stdev' in checkbox_values and 'show-pfam' in checkbox_values:
            protein_fig = protein_graphic_plotly(protein_list_cutoff, difference_metric='area_sum', show_pfam=True, show_stdev = True)
        elif checkbox_values and 'show-stdev' in checkbox_values:
            protein_fig = protein_graphic_plotly(protein_list_cutoff, difference_metric='area_sum', show_stdev = True)
        elif checkbox_values and 'show-pfam' in checkbox_values:
            protein_fig = protein_graphic_plotly(protein_list_cutoff, difference_metric='area_sum', show_pfam = True)
        else:
            protein_fig = protein_graphic_plotly(protein_list_cutoff, difference_metric='area_sum')
        return protein_fig, triv_names, datatable
    
    else:
        return {}, [], html.Div()

peptide_lists=[]
def create_peptide_fig(clickData, search_protein, n_clicks_sum, n_clicks_mean):
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
        if n_clicks_sum == 0 and n_clicks_mean == 0 or n_clicks_sum > n_clicks_mean:
            peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=False)
        else:
            peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=True)
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
                'minWdith':95,
            },
            style_table={'height': '400px', 'width':'500px', 'overflowY': 'auto', 'overflowX':'auto'}
    )   
        return peptide_fig, search_text, html.Div(datatable)
    
    else:
        return {}, search_text, []


def amino_acid_dropdown(n_clicks_complete_proteome, n_clicks_selected_protein, radioitem_value):
    if n_clicks_complete_proteome > n_clicks_selected_protein and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(protein_lists[-1], peptide_or_protein_list = 'protein_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Complete proteome'
    if n_clicks_selected_protein > n_clicks_complete_proteome and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_lists[-1], peptide_or_protein_list = 'peptide_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Selected Protein'
    else:
        return  {},{},{}, {},{},{}, ''


def generate_hover_graphs(hoverData):
    if hoverData:
        accession = hoverData['points'][0]['customdata'][-1]
        protein_list = protein_lists[-1]
        fig = all_sample_bar_chart(protein_list, accession=accession, metric='area_sum',)
        return fig
    else:
        return {}


app.callback(
    Output('hover-all-protein-samples', 'figure'),
    Input('protein-fig','hoverData')
)(generate_hover_graphs)

app.callback(
    Output('complete-aa-seq-fig-g1', 'figure'),
    Output('first-aa-fig-g1', 'figure'),
    Output('last-aa-fig-g1', 'figure'),
    Output('complete-aa-seq-fig-g2', 'figure'),
    Output('first-aa-fig-g2', 'figure'),
    Output('last-aa-fig-g2', 'figure'),
    Output('complete-or-selected', 'children'),
    Input('complete-proteome', 'n_clicks_timestamp'),
    Input('selected-protein','n_clicks_timestamp'),
    Input('aa-radioitems', 'value')
)(amino_acid_dropdown)

app.callback(
    Output('peptide-fig', 'figure'),
    Output('search-protein', 'value'),
    Output('peptide-info-table', 'children'),
    Input('protein-fig', 'clickData'),
    Input('search-protein', 'value'),
    Input('peptide-dropdown-sum', 'n_clicks_timestamp'),
    Input('peptide-dropdown-mean','n_clicks_timestamp'),
)(create_peptide_fig)

app.callback(
    Output("output-filename-1", "data"),
    Output("sample-collapse-1", "data"),
    [Input("upload-data-1", "contents"), Input("upload-data-1", "filename")],
)(update_file_list)

app.callback(
    Output("output-filename-2", "data"),
    Output("sample-collapse-2", "data"),
    [Input("upload-data-2", "contents"), Input("upload-data-2", "filename")],
)(update_file_list)

app.callback(
    Output('protein-fig', 'figure'),
    Output('protein-list', 'children'),
    Output('protein-info-table','children'),
    Input("close-modal-file", "n_clicks_timestamp"),
    Input('protein-checklist','value'),
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
    Output("how-to-use-collapse", "is_open"),
    [Input("how-to-use-collapse-button", "n_clicks")],
    [State("how-to-use-collapse", "is_open")],
)(toggle_collapse)

app.callback(
    Output("sample-collapse", "is_open"),
    [Input("sample-collapse-button", "n_clicks")],
    [State("sample-collapse", "is_open")],
)(toggle_collapse)




if __name__ == '__main__':
    app.run_server(debug=True)
