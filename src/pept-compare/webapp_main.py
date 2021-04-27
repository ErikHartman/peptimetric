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
from dash import no_update
import plotly.graph_objects as go

from methods import protein_list_to_json, json_to_protein_list, peptide_list_to_json, json_to_peptide_list
from methods import make_peptide_dfs, concatenate_dataframes, merge_dataframes, create_protein_list, create_protein_df_fig, create_protein_fig , create_peptide_list, stacked_samples_peptide
from methods import amino_acid_piecharts, all_sample_bar_chart, create_peptide_list_from_trivname
from methods import apply_protein_cutoffs, apply_peptide_cutoffs, get_unique_and_common_proteins, create_venn_bar
from methods import proteins_present_in_all_samples, create_protein_datatable, create_peptide_datatable, log_intensity, normalize_data, create_length_histogram
from texts_for_webapp import how_to_use, Documentation, contact_text


app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(id='page-content', fluid=True, className='vh-100'),
 ])

#---------------------------------------PAGE-ELEMENTS------------------------------------------------
file_columns = ['Sample', 'File']

modal_file = html.Div([

    dbc.Button("Upload files", id="open-modal-file", color='secondary',  outline=True, style={'border-color':'transparent'}, className="mr-1"),
        dbc.Modal([
                dbc.ModalHeader("Upload files", className="font-weight-bold"),
                    dbc.Row([
                        dbc.Col(dbc.ModalBody('Group 1', className='ml-auto text-center font-weight-bold')),
                        dbc.Col(dbc.ModalBody('Group 2', className='ml-auto text-center font-weight-bold')),
                        
                    ]),
                    dbc.Row([
                        dbc.Col(dcc.Upload(id='upload-data-1', children=dbc.Button('Select files', style={'padding':10}), multiple=True), className="text-center ml-auto"),
                        dbc.Col(dcc.Upload(id='upload-data-2', children=dbc.Button('Select files', style={'padding':10}), multiple=True), className="text-center ml-auto"),
                        
                    ]),    
                    dbc.Row([
                        dbc.Col(dash_table.DataTable(
                                id = 'output-filename-1',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                        'if' : {'row_index':'odd'},
                                        'backgroundColor' : '#92dbb0'
                                    }
                                    ],
                                    style_header={
                                        'textAlign':'center',
                                        'fontWeight': 'bold',
                                        'font-family':'Roboto'
                                    },
                                    style_cell={
                                        'textAlign':'left',
                                        'padding':'5px',
                                        'font-family':'Roboto',
                                        'fontSize':12,
                                    },
                                    style_table = {'padding':10}
                            )),
                        dbc.Col(dash_table.DataTable(
                                id = 'output-filename-2',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                        'if' : {'row_index':'odd'},
                                        'backgroundColor' : '#92dbb0'
                                    }
                                    ],
                                    style_header={
                                        'textAlign':'center',
                                        'fontWeight': 'bold',
                                        'font-family':'Roboto'
                                    },
                                    style_cell={
                                        'textAlign':'left',
                                        'padding':'5px',
                                        'font-family':'Roboto',
                                        'fontSize':12,
                                    },
                                    style_table = {'padding':10}
                            )),
                    ]),
                dbc.ModalFooter([
                    dbc.Label('*upload 3 or more files for statistical analysis'),
                    dbc.Button("Upload files", color = 'primary', id="close-modal-file", className="ml-auto", n_clicks_timestamp=0)
                ]),
            ],
            id="modal-file",
            centered=True,
            scrollable=True,
        )])


protein_tab = dbc.Form([
                    dbc.Col([dbc.Label("Intensity", style={'padding':10, 'margin':5}),
                                    dbc.Label('Spectral count', style={'padding':10}),
                                    dbc.Label('Number of peptides', className='mr-2', style={'padding':10})]),
                        dbc.Col([dbc.Input(placeholder='0', type='number', className='ml-auto', min=0, step=0.00000001, 
                        id='tot-intensity-cutoff', value=0, style={'padding':10, 'margin':5}),
                        dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, step=0.00000001, 
                        id='tot-spc-cutoff', value=0, style={'padding':10, 'margin-bottom':5}),
                        dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, 
                        id='nbr-of-peptides-cutoff', value=0, style={'padding':10})]),
                    dbc.FormGroup([
                        dbc.Col(dbc.Checklist(
                         options=[
                            {"label": "Only show proteins present in all samples", "value": 'present-in-all-samples'},
                            ],
                            value=[],
                            id='proteins-present-in-all-samples-checkbox'), width = 12)],
                     className="mr-3") 
                ],
                inline=True),

peptide_tab = dbc.Form([
                        dbc.Col([dbc.Label("Intensity", style={'padding':10, 'margin':5}),
                                    dbc.Label('Spectral count', style={'padding':10})]),
                        dbc.Col([dbc.Input(placeholder='0', type='number', className='ml-auto', min=0, step=0.001, 
                        id='peptide-intensity-cutoff', value=0, style={'padding':10, 'margin':5}),
                        dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, step=0.001, 
                        id='peptide-spc-cutoff', value=0, style={'padding':10})]),
                    dbc.FormGroup([
                        dbc.Col(dbc.Checklist(
                            options=[
                                {"label": "Remove RT outliers", "value": 'RT'},
                                {"label": "Remove CCS outliers", "value": 'CCS'},
                                ],
                                value=[], style={'padding':10},
                                id='RT-CCS-checkbox'), width=12)],
                     className="mr-3")
                ],
                inline=True),

modal_cutoff = dbc.Modal([
                dbc.ModalHeader("Cutoff settings", className="font-weight-bold"),
                dbc.Tabs([
                    dbc.Tab(peptide_tab, label='Peptide'),
                    dbc.Tab(protein_tab, label='Protein'),
                    
                ]),
                dbc.ModalFooter([
                    dbc.Label('*cutoffs will be applied after normalization'),
                    dbc.Button("Apply", id="close-modal-cutoff", color="primary", className="ml-auto")
                ]),
            ],
            id="modal-cutoff",
            size='m',
            centered=True,
              )



modal_feedback = html.Div([
    dbc.Button("Feedback", id="open-modal-feedback", color='secondary', outline=True, style={'border-color':'transparent'}, className="mr-1"),
    dbc.Modal([
        dbc.ModalHeader('Feedback'),
                dbc.ModalBody(
                    html.Div(
                        html.P(contact_text)
                    )
                ),
                
                dbc.ModalFooter([
                    dbc.Button("Close", id="close-modal-feedback", color='primary', className="ml-auto")
                ]),
            ],
            id="modal-feedback",
            size='m',
            centered=True,
              )])


normalization_modal = dbc.Modal([
                dbc.ModalHeader("Normalize data", className="font-weight-bold", style={'padding':10}),
                dbc.ModalBody([
                    dbc.FormGroup([
                        dbc.Col(dbc.RadioItems(
                                    options=[
                                    {'label': 'Normalize on global intensity', 'value': 'global-intensity'},
                                    {'label': 'Normalize on housekeeping protein', 'value': 'housekeeping-protein'}
                                    ],
                                    value='',
                                    id='normalization-radioitems',
                                    style={'padding':10}
                        )),
                        dbc.Col(dbc.Input(
                                    id='housekeeping-protein-input',
                                    placeholder='Protein name (e.g. ALBU_HUMAN)',
                                    name = 'text',
                                    debounce=True,
                                    inputMode='latin',
                                    minLength=0, maxLength=30,
                                    size = '10',
                                    list = 'protein-list',
                                    className="ml-auto",
                                    disabled=True
                                ), width=8,
                        )
                    ]),
                dbc.Checklist(
                    options=[
                        {"label": "I've already taken the log of my intensities", "value": True},
                    ],
                    value=[],
                    id="log-checkbox",
                ),
                ]),
                
                dbc.ModalFooter(
                    dbc.Button("Apply", id="close-modal-normalization", color="primary", className="ml-auto")
                ),
            ],
            id="modal-normalization",
            size='m',
            centered=True,
              )


navbar = dbc.Navbar(
    [
        dbc.NavLink("Peptimetric", href = '/', style = {'color':'grey', 'font-size':20,  'font-weight':'bold', 'font':'Roboto', 'text-transform':'lowercase'} ),
        modal_file,
        dbc.Button('Normalization', id="open-modal-normalization", color='secondary', outline=True, style={'border-color':'transparent'}, className='mr-1'),
        normalization_modal,
        dbc.Button("Cutoffs", id="open-modal-cutoff", color='secondary', outline=True, style={'border-color':'transparent'}, className='mr-1'),
        modal_cutoff,
        dbc.Nav([
        modal_feedback,
        dbc.NavLink('Documentation', href='/documentation'),
        dbc.NavLink('Home', href='/')
        ],
        navbar=True,
        className="ml-auto",)
    ],   
)

amino_acid_pie_dropdown = dcc.Dropdown(
    id= 'amino-acid-pie-dropdown',
    placeholder='Select view',
    value='',
    options=[
        {'label': 'Selected protein', 'value': 'selected-protein'},
        {'label': 'Complete proteome', 'value': 'complete-proteome'},
    ],
    
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
            disabled=True,
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
            dbc.CardBody(
                how_to_use, style={"maxHeight": "300px", "overflowY": "scroll"}
            ),
            dbc.CardFooter(dbc.Button('Load example data', id='load-example-data', color='primary')
            ),
            
            ]),
            id="how-to-use-collapse",
            className="mb-4",
            is_open=True,
            
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
                dbc.Col('Group 1', className='font-weight-bold text-center'),
                dbc.Col('Group 2', className='font-weight-bold text-center')
                ]),
            dbc.Row([
            dbc.Col(dash_table.DataTable(
                                id = 'sample-collapse-1',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                        'if' : {'row_index':'odd'},
                                        'backgroundColor' : '#92dbb0'
                                    }
                                    ],
                                    style_header={
                                        'textAlign':'center',
                                        'fontWeight': 'bold',
                                        'font-family':'Roboto'
                                    },
                                    style_cell={
                                        'textAlign':'left',
                                        'padding':'5px',
                                        'font-family':'Roboto',
                                        'fontSize':12,
                                    },
                                    style_table = {'padding':10} 
                            )),
            dbc.Col(dash_table.DataTable(
                                id = 'sample-collapse-2',
                                columns=[{"name": i, "id": i} for i in file_columns],
                                data=[],
                                style_data_conditional = [{
                                        'if' : {'row_index':'odd'},
                                        'backgroundColor' : '#92dbb0'
                                    }
                                    ],
                                    style_header={
                                        'textAlign':'center',
                                        'fontWeight': 'bold',
                                        'font-family':'Roboto'
                                    },
                                    style_cell={
                                        'textAlign':'left',
                                        'padding':'5px',
                                        'font-family':'Roboto',
                                        'fontSize':12,
                                    },
                                    style_table = {'padding':10} )) 
                ]),
            ]),
            id="sample-collapse",
            className="mb-4"
        ),
    ]
)

protein_fig_radioitems = html.Div([
    dcc.Dropdown(
        placeholder = 'Select difference metric...',
        options=[
        {'label': 'Intensity sum', 'value': 'area_sum'},
        {'label': 'Intensity mean', 'value': 'area_mean'},
        {'label': 'SpC sum', 'value': 'spc_sum'},
        {'label': 'SpC mean', 'value': 'spc_mean'}
        ],
        value='area_sum',
        id='protein-radioitems',
    )
])

protein_fig = html.Div([
    dbc.Row([
        html.Img(src = app.get_asset_url ('scatter.png'), style={'height':'6%', 'width':'6%'}),
        html.H3('Protein View'),
    ]),
        dbc.Row([
            dbc.Col(search_protein, width={'size':3}),
            dbc.Col(dbc.Checklist(
                className='ml-auto',
                id='protein-checklist',
                switch=True,
                inline=True,
                options=[
                    {'label': 'Show standard deviation', 'value': 'show-stdev'},

                    ],
                )                 
            ),
            dbc.Col(protein_fig_radioitems),
            dbc.Col(dbc.Button('Generate protein graph', id='generate-protein-graph', color='success'))
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='protein-fig', figure={}, config={'toImageButtonOptions':{'format': 'jpeg'}, 'displaylogo':False})
        )
        
        ])

all_samples_protein_fig = html.Div([
    dcc.Graph(id='hover-all-protein-samples', figure={}, style={'height': 300, 'width':500},
    config={'displayModeBar': False}
)])

peptide_fig_radioitems = html.Div([
    dbc.RadioItems(
        options=[
        {'label': 'Intensity', 'value': 'area'},
        {'label': 'Spectral Count', 'value': 'spectral_count'}
        ],
        value='area',
        id='peptide-radioitems',
        
    )
])

peptide_fig_radioitems_sum_or_mean = html.Div([ 
    dbc.RadioItems(
        options=[
        {'label': 'View all samples', 'value': False},
        {'label': 'View mean of all samples', 'value': True}
        ],
        value=False,
        id='sum-or-mean-radio',
       
    )
])

peptide_fig = html.Div([
        dbc.Row([
        html.Img(src = app.get_asset_url ('bar.png'), style={'height':'6%', 'width':'6%'}),
        html.H3('Peptide View'),
    ]),
        dbc.Row([
            dbc.Col(peptide_fig_radioitems_sum_or_mean),
            dbc.Col(peptide_fig_radioitems)    ,
            dbc.Col(dbc.Button('Choose a protein', disabled=True, id='generate-peptide-fig', color='success'))             
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='peptide-fig', figure={}, config={'displaylogo': False})
        ),
        ])

amino_acid_radioitems = html.Div([
        dbc.RadioItems(
            options=[
                {"label": "Intensity", "value": 'area'},
                {"label": "Spectral Count", "value": 'spectral_count'}, 
            ], 
            value='area',
            id="aa-radioitems", 
            inline=True,

        )
])

amino_acid_figs = html.Div([
        html.H3('Amino Acid Profile'),
        dbc.Row(dbc.Col(amino_acid_radioitems)),
        dcc.Loading(type='cube', color = '#76b382',
            children=[ dbc.Row([
                dbc.Col(dcc.Graph(id='aa-fig', figure={}, style={'height': '700px'}, config={'displayModeBar': False})),
            ]) 
            ]
        )
    ])

peptide_length_fig = html.Div([
    html.H3('Peptide Length'),
    dcc.Loading(type='cube', color = '#76b382',
        children=[    
        dbc.Row([
            dbc.Col(dcc.Graph(id='peptide-length-fig', figure={}, config={'displayModeBar': False}))
        ]),])
])

venn_bar_fig = html.Div([
    html.H3('Peptide overlap'),
    dcc.Loading(type='cube', color = '#76b382',
        children=[    
        dbc.Row([
            dbc.Col([
                dcc.Graph(id='venn-bar', figure={}, config={'displayModeBar': False})]),
        ])]
    )
])

start_table_df = pd.DataFrame(columns=['No data'])
protein_info = html.Div(dash_table.DataTable(
            data=start_table_df.to_dict('records'), 
            columns = [{'id': c, 'name': c} for c in start_table_df.columns],
            id='protein-info-table',
            sort_action='native',
            fixed_rows={'headers': True},
            filter_action='native',
            virtualization=True,
            row_selectable="multi",
            export_format='xlsx',
            selected_rows=[],
            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : '#92dbb0'
            }
            ],
            style_header={
                'textAlign':'center',
                'fontWeight': 'bold',
                'font-family':'Roboto'
            },
            style_cell={
                'textAlign':'left',
                'padding':'5px',
                'maxWidth': 105,
                'minWidth': 105,
                'font-family':'Roboto',
                'fontSize':12,
            },
            style_table={'height': '200px', 'width':'500px', 'overflowY': 'auto','overflowX':'auto'}
    ),
)


peptide_info = html.Div(dash_table.DataTable(id='peptide-info-table',
            data=start_table_df.to_dict('records'), 
            columns = [{'id': c, 'name': c} for c in start_table_df.columns],
            sort_action='native',
            fixed_rows={'headers': True},
            filter_action='native',
            virtualization=True,
            export_format='xlsx',

            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : '#92dbb0'
            }
            ],
            style_header={
                'textAlign':'center',
                'fontWeight': 'bold',
                'font-family':'Roboto'
            },
            style_cell={
                'textAlign':'left',
                'padding':'5px',
                'maxWidth': 105,
                'minWidth': 105,
                'font-family':'Roboto',
                'fontSize':12,
            },
            style_table={'height': '400px', 'width':'500px', 'overflowY': 'auto','overflowX':'auto'}
    ),
)

documentation = dbc.Col([
    dbc.Row([html.Img(src = app.get_asset_url ('document.png'), style={'height':'4%', 'width':'4%'}),
        dbc.Col(html.H1('Documentation'))]),
    html.Hr(),
    Documentation,
    
])


hidden_divs = html.Div([
    dcc.Store(id='cutoff-value-holder'),
    dcc.Store(id='protein-list-df-holder'),
    dcc.Store(id='peptide-list-df-holder'),
    dcc.Store(id='df_g1-holder'),
    dcc.Store(id='df_g2-holder'),
    dcc.Store(id='peptide-data-holder'),
    dcc.Store(id='protein-fig-holder'),
    dcc.Store(id='normalization-holder'),
    dcc.Loading(color = '#76b382', style={'backgroundColor': 'transparent'}, className = 'loader-wrapper', fullscreen = True, type = 'default', id='process-data-loading', children = [dcc.Store(id='protein-datatable-holder')]),
    dcc.Store(id='housekeeping-protein-holder'),
    
])
#---------------------------PAGES---------------------------------------------------------------
main_page = dbc.Container([
    
    dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4", style={'padding-left':0, 'padding-right':0},)
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
        html.Img(src = app.get_asset_url ('pie.png'), style={'height':'4%', 'width':'4%'}),
        html.H3('General characteristics', style={'bottom-margin':0}),
    ]),
    dbc.Row([dbc.Col(amino_acid_pie_dropdown, width=2, style={'padding':15})]),
    dbc.Row([
        dbc.Col(peptide_length_fig, width={'size':8}),
        dbc.Col(venn_bar_fig, width={'size':4}),
    ]),  
    dbc.Row(dbc.Col(amino_acid_figs)),
    hidden_divs,
], fluid=True, style={'padding':0})

documentation_page = dbc.Container([
     dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4", style={'padding-right':0, 'padding-left':0})
    ]),
    dbc.Row([
        documentation,
    ])

    
], fluid=True, style={'padding':0})

#-----------------DEFS AND CALLBACKS--------------------------------------------------------------
def display_page(pathname):
    if pathname == '/':
        return main_page
    elif pathname == '/documentation':
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


def update_file_list(contents, filename, log_checkbox):
    if filename:
        file_list = []
        i=0
        for f in filename:
            s = [f"S{i}", f]
            file_list.append(s)
            i+=1  
        master_df = update_data_frame(contents, filename, log_checkbox)
        df = pd.DataFrame(file_list, columns = ['Sample', 'File'])
        return df.to_dict('rows'), df.to_dict('rows'), master_df.to_json()
    else:
        return [], [], []

def update_data_frame(contents, filename, log_checkbox):
    decoded_list = []
    for f, content in zip(filename, contents):
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        decoded_list.append(io.BytesIO(decoded))
    dfs = make_peptide_dfs(decoded_list, filename)
    master_df = concatenate_dataframes(dfs)
    if log_checkbox:
        return master_df
    master_df = log_intensity(master_df)
    return master_df


def set_cutoffs(tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT_CCS_checkbox, proteins_present_in_all_samples):
    RT = False
    CCS = False
    present_in_all_samples = False
    if 'RT' in RT_CCS_checkbox:
        RT=True
    if 'CCS' in RT_CCS_checkbox:
        CCS=True
    if 'present-in-all-samples' in proteins_present_in_all_samples:
        proteins_present_in_all_samples = True
    return [tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CCS, present_in_all_samples]


def process_protein_data(apply_normalization_n_clicks, n_clicks_close_file, apply_cutoffs_button, cutoff_values, df_g1, df_g2, radioitems_normalization, housekeeping_protein):
    if apply_cutoffs_button:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CCS, present_in_all_samples = cutoff_values
    else:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CCS, present_in_all_samples = 0,0,0,0,0,False,False, False
    triv_names = []
    protein_list = []
    protein_list_cutoff = []
    if n_clicks_close_file and df_g1 and df_g2:
        g1 = pd.read_json(df_g1)
        g2 = pd.read_json(df_g2)
        master = merge_dataframes(g1,g2)
        protein_list = create_protein_list(master)
        if 'global-intensity' in radioitems_normalization:
            protein_list = normalize_data(protein_list, housekeeping_protein=False)
        elif 'housekeeping-protein' in radioitems_normalization and housekeeping_protein != '':
            protein_list = normalize_data(protein_list, housekeeping_protein = housekeeping_protein)
        protein_list_cutoff = apply_peptide_cutoffs(protein_list, area=pep_intensity_co, spc=pep_spc_co, rt=RT, ccs=CCS)
        protein_list_cutoff = apply_protein_cutoffs(protein_list_cutoff, nbr_of_peptides=nbr_of_peptides_co, tot_area=tot_intensity_co, tot_spc=tot_spc_co)
        if len(protein_list_cutoff) < 1:
            return [], [], [], []
        if present_in_all_samples:
            protein_list_cutoff = proteins_present_in_all_samples(protein_list_cutoff)
        
        protein_list_json = protein_list_to_json(protein_list_cutoff)
        df_fig = create_protein_df_fig(protein_list_cutoff)
        df_protein_info = create_protein_datatable(protein_list_cutoff)
        df_protein_info.fillna(0, inplace=True)
        if len(protein_list) > 1:
            for protein in protein_list_cutoff:
                triv_names.append(html.Option(value=protein.get_trivial_name()))
        return triv_names, df_fig.to_json(), df_protein_info.to_json(), protein_list_json

    else:
        return [], [], [], []

def create_protein_figure_and_table(rows, derived_virtual_selected_rows, search_protein, clickData, protein_radioitems_value, checkbox_values, generate_protein_graph_n_clicks, df_fig, df_protein_info, protein_fig, protein_list_json):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    highlighted_triv_names = ['Choose protein']
    disabled = True
    if df_protein_info:
        df_protein_info = pd.read_json(df_protein_info)
    if clickData or search_protein or derived_virtual_selected_rows:
        disabled = False
        protein_fig = go.Figure(protein_fig)
        highlighted_triv_names = []
        triv_names_holder = df_protein_info['Protein'].array
        if str(changed_id) == 'search-protein.value' and search_protein in triv_names_holder:
            highlighted_triv_names.append(search_protein)
        elif str(changed_id) == 'protein-fig.clickData':
            highlighted_triv_names.append(clickData['points'][0]['customdata'][0])
        elif str(changed_id) == 'protein-info-table.derived_virtual_data' or str(changed_id) == 'protein-info-table.derived_virtual_selected_rows':
            selected_rows_df = pd.DataFrame(rows)
            highlighted_triv_names = list(selected_rows_df.iloc[derived_virtual_selected_rows, 0])
            if not rows:
                highlighted_triv_names = []

        marker_color_list = ['rgba(0,0,0,0)' for n in range(len(triv_names_holder))]
        for triv_name in highlighted_triv_names:
            for i in range(len(triv_names_holder)):
                if triv_name == str(triv_names_holder[i]):
                    marker_color_list[i] = 'rgba(242, 89, 0, 1)'
        protein_fig.update_traces(marker=dict(line=dict(width=3, color=marker_color_list)),selector=dict(mode='markers'))
    if df_fig:
        df_fig = pd.read_json(df_fig)
        protein_list = json_to_protein_list(protein_list_json)
        if 'area_sum' in protein_radioitems_value:
            difference_metric = 'area_sum'
            columns = ['Protein','UniProt ID','#peptides_g1','#peptides_g2','intensity_sum_g1', 'intensity_sum_g1_sd','intensity_sum_g2', 'intensity_sum_g2_sd', 'intensity_p_value_sum']
            sort = ['intensity G2', 'intensity G1']
            rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2', 'intensity_sum_g1': 'intensity G1', 'intensity_sum_g2':'intensity G2', 'intensity_sum_g1_sd':'SD', 'intensity_sum_g2_sd': 'SD', 'intensity_p_value_sum':'p-value'
            }
        elif 'area_mean' in protein_radioitems_value:
            difference_metric = 'area_mean'
            columns = ['Protein','UniProt ID','#peptides_g1','#peptides_g2','intensity_mean_g1', 'intensity_mean_g1_sd','intensity_mean_g2', 'intensity_mean_g2_sd', 'intensity_p_value_mean']
            sort = ['intensity G2', 'intensity G1']
            rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','intensity_mean_g1': 'intensity G1', 'intensity_mean_g2':'intensity G2', 'intensity_mean_g1_sd':'SD', 'intensity_mean_g2_sd': 'SD', 'intensity_p_value_mean':'p-value'
            }
        elif 'spc_sum' in protein_radioitems_value:
            difference_metric = 'spc_sum'
            columns = ['Protein','UniProt ID','#peptides_g1','#peptides_g2','spc_sum_g1', 'spc_sum_g1_sd','spc_sum_g2', 'spc_sum_g2_sd', 'spc_p_value_sum']
            sort = ['SpC G2','SpC G1']
            rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','spc_sum_g1': 'SpC G1', 'spc_sum_g2':'SpC G2', 'spc_sum_g1_sd':'SD', 'spc_sum_g2_sd': 'SD', 'spc_p_value_sum':'p-value'
            }
        else:
            difference_metric = 'spc_mean'
            columns = ['Protein','UniProt ID','#peptides_g1','#peptides_g2','spc_mean_g1', 'spc_mean_g1_sd','spc_mean_g2', 'spc_mean_g2_sd', 'spc_p_value_mean']
            sort = ['SpC G2','SpC G1']
            rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','spc_mean_g1': 'SpC G1', 'spc_mean_g2':'SpC G2', 'spc_mean_g1_sd':'SD', 'spc_mean_g2_sd': 'SD', 'spc_p_value_mean':'p-value'
            }
        
        df_protein_info = df_protein_info[columns]
        df_protein_info.rename(columns = rename, inplace=True)
        df_protein_info.sort_values(by=sort, ascending=False, inplace=True)
        protein_info_data = df_protein_info.to_dict('rows')
        protein_info_columns=[{"name": str(i), "id": str(i)} for i in df_protein_info.columns]
        if str(changed_id) == 'generate-protein-graph.n_clicks' or str(changed_id) == 'protein-radioitems.value' or str(changed_id) == 'protein-checkbox.value':
            if checkbox_values and 'show-stdev' in checkbox_values :
                protein_fig = create_protein_fig(df_fig, protein_list, show_stdev = True,  difference_metric = difference_metric)
            else:
                protein_fig = create_protein_fig(df_fig, protein_list, difference_metric = difference_metric)
        if len(highlighted_triv_names) < 1:
            highlighted_triv_names = ['Choose protein']
            disabled = True
        return protein_fig, protein_info_data, protein_info_columns, disabled, str(highlighted_triv_names[0])
    else:
        return {}, start_table_df.to_dict('records'), [{'id': '', 'name': ''}], disabled, ['Choose protein']

@app.callback(
    Output('normalization-holder', 'children'),
    Output('housekeeping-protein-holder', 'children'),
    Input('normalization-radioitems', 'value'),
    Input('housekeeping-protein-input','value'),
)
def get_normalization_data(radioitems_normalization, housekeeping_protein):
    return radioitems_normalization, housekeeping_protein


def create_peptide_fig(n_clicks_generate_peptide_fig, sum_or_mean_radio, peptide_radioitems_value, button_label, protein_list_json, peptide_list_json):
    if peptide_radioitems_value == 'area':
        columns = ['Peptide','Start','End','intensity_g1','intensity_g1_sd','intensity_g2', 'intensity_g2_sd']
        sort = ['intensity G1','intensity G2']
        rename = {'intensity_g1':'intensity G1', 'intensity_g1_sd':'SD','intensity_g2':'intensity G2','intensity_g2_sd':'SD'}
    else:
        columns = ['Peptide','Start','End','spc_g1', 'spc_g1_sd','spc_g2', 'spc_g2_sd']
        sort = ['SpC G1','SpC G2']
        rename= {'spc_g1':'SpC G1', 'spc_g1_sd':'SD','spc_g2':'SpC G2','spc_g2_sd':'SD'}


    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if changed_id == 'generate-peptide-fig.n_clicks':
        protein_list = json_to_protein_list(protein_list_json)
        trivname = button_label.split(' ')[-1]
        peptide_list = create_peptide_list_from_trivname(protein_list, str(trivname))
        peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=sum_or_mean_radio, difference_metric=peptide_radioitems_value)
        df_peptide_info = create_peptide_datatable(peptide_list)
        df_peptide_info.fillna(0, inplace=True)
        df_peptide_info = df_peptide_info[columns]
        df_peptide_info.rename(columns=rename, inplace=True)
        df_peptide_info.sort_values(by=sort, ascending=False, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows')
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in df_peptide_info.columns]
        return peptide_fig,  peptide_table_data, peptide_table_columns, peptide_list_to_json(peptide_list)
    elif peptide_list_json and (changed_id == 'sum-or-mean-radio.value' or 'peptide-radioitems.value'):
        peptide_list = json_to_peptide_list(peptide_list_json)
        peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=sum_or_mean_radio, difference_metric=peptide_radioitems_value)
        df_peptide_info = create_peptide_datatable(peptide_list)
        df_peptide_info.fillna(0, inplace=True)
        df_peptide_info = df_peptide_info[columns]
        df_peptide_info.rename(columns=rename, inplace=True)
        df_peptide_info.sort_values(by=sort, ascending=False, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows')
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in df_peptide_info.columns]
        return peptide_fig, peptide_table_data, peptide_table_columns, peptide_list_json
    
    else:
        return {}, start_table_df.to_dict('records'), [{'id': '', 'name': ''}], []


def amino_acid_dropdown(dropdown_values, radioitem_value, protein_list_json, peptide_list_json):
    if dropdown_values and 'complete-proteome' in dropdown_values and protein_list_json:
        protein_list = json_to_protein_list(protein_list_json)
        fig = amino_acid_piecharts(protein_list, peptide_or_protein_list = 'protein_list', difference_metric = radioitem_value)
        return fig
    elif dropdown_values and 'selected-protein' in dropdown_values and peptide_list_json:
        peptide_list = json_to_peptide_list(peptide_list_json)
        fig = amino_acid_piecharts(peptide_list, peptide_or_protein_list = 'peptide_list', difference_metric = radioitem_value)
        return fig
    else:
        return  {}


def generate_hover_graphs(hoverData, protein_radioitems_value, protein_list_json):
    if hoverData:
        accession = hoverData['points'][0]['customdata'][-1]
        protein_list = json_to_protein_list(protein_list_json)
        fig = all_sample_bar_chart(protein_list, accession=accession, metric=protein_radioitems_value)
        return fig
    else:
        return {}

def enable_input_housekeeping_protein(normalization_radio):
    if 'housekeeping-protein' in normalization_radio:
        return False
    else:
        return True

def enable_input_search_protein(protein_list):
        if not protein_list:
            return True
        else:
            return False

def enable_generate_protein_graph(protein_list):
        if not protein_list:
            return True
        else:
            return False

def create_peptide_length_dropdown(length_dropdown_values, protein_list_json, peptide_list_json):
    
    if length_dropdown_values and 'complete-proteome'  in length_dropdown_values and protein_list_json:
        protein_list = json_to_protein_list(protein_list_json)
        length_fig = create_length_histogram(protein_list, peptide_or_protein_list='protein_list')
        return length_fig
    elif length_dropdown_values and 'selected-protein' in length_dropdown_values and peptide_list_json:
        peptide_list = json_to_peptide_list(peptide_list_json)
        length_fig= create_length_histogram(peptide_list, peptide_or_protein_list='peptide_list')
        return length_fig
    else:
        return {}

def create_venn_bar_fig(length_dropdown_values, protein_list_json, peptide_list_json):
    if length_dropdown_values and 'complete-proteome'  in length_dropdown_values and protein_list_json:
        protein_list = json_to_protein_list(protein_list_json)
        venn_bar = create_venn_bar(protein_list, complete_proteome=True)
        return venn_bar
    elif length_dropdown_values and 'selected-protein' in length_dropdown_values and peptide_list_json:
        peptide_list = json_to_peptide_list(peptide_list_json)
        venn_bar= create_venn_bar(peptide_list, complete_proteome=False)
        return venn_bar
    else:
        return {}



app.callback(
    Output('generate-protein-graph', 'disabled'),
    Input('protein-list', 'children')
)(enable_generate_protein_graph)

app.callback(
    Output('search-protein', 'disabled'),
    Input('protein-list', 'children')
)(enable_input_search_protein)
    

app.callback(
    Output('housekeeping-protein-input', 'disabled'),
    Input('normalization-radioitems', 'value')
)(enable_input_housekeeping_protein)
    


app.callback(
    Output('hover-all-protein-samples', 'figure'),
    Input('protein-fig','hoverData'),
    Input('protein-radioitems', 'value'),
    State('protein-list-df-holder','children'),
)(generate_hover_graphs)

app.callback(
    Output('aa-fig', 'figure'),

    Input('amino-acid-pie-dropdown', 'value'),
    Input('aa-radioitems', 'value'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(amino_acid_dropdown)

app.callback(
    Output('peptide-length-fig', 'figure'),
    Input('amino-acid-pie-dropdown', 'value'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(create_peptide_length_dropdown)

app.callback(
    Output('venn-bar', 'figure'),
    Input('amino-acid-pie-dropdown', 'value'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(create_venn_bar_fig)


app.callback(
    Output('peptide-fig', 'figure'),
    Output('peptide-info-table', 'data'),
    Output('peptide-info-table', 'columns'),
    Output('peptide-list-df-holder', 'children'),
    Input('generate-peptide-fig', 'n_clicks'),
    Input('sum-or-mean-radio', 'value'),
    Input('peptide-radioitems', 'value'),
    State('generate-peptide-fig', 'children'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(create_peptide_fig)

app.callback(
    Output("output-filename-1", "data"),
    Output("sample-collapse-1", "data"),
    Output('df_g1-holder', 'children'),
    [Input("upload-data-1", "contents"), Input("upload-data-1", "filename"),
    Input('log-checkbox', 'value')],
)(update_file_list)

app.callback(
    Output("output-filename-2", "data"),
    Output("sample-collapse-2", "data"),
    Output('df_g2-holder','children'),
    [Input("upload-data-2", "contents"), Input("upload-data-2", "filename"),
    Input('log-checkbox', 'value')],
)(update_file_list)

app.callback(
    Output('protein-list', 'children'),
    Output('protein-fig-holder', 'children'),
    Output('protein-datatable-holder','children'),
    Output('protein-list-df-holder', 'children'),
    Input('close-modal-normalization', 'n_clicks'),
    Input("close-modal-file", "n_clicks_timestamp"),
    Input('close-modal-cutoff', 'n_clicks'),
    State('cutoff-value-holder', 'children'),
    State('df_g1-holder', 'children'),
    State('df_g2-holder', 'children'),
    State('normalization-holder', 'children'),
    State('housekeeping-protein-holder', 'children'),
    )(process_protein_data)


app.callback(
    Output('protein-fig', 'figure'),
    Output('protein-info-table', 'data'),
    Output('protein-info-table', 'columns'),
    Output('generate-peptide-fig', 'disabled'),
    Output('generate-peptide-fig', 'children'),
    Input('protein-info-table', "derived_virtual_data"),
    Input('protein-info-table', 'derived_virtual_selected_rows'),
    Input('search-protein', 'value'),
    Input('protein-fig', 'clickData'),
    Input('protein-radioitems','value'),
    Input('protein-checklist', 'value'),
    Input('generate-protein-graph', 'n_clicks'),
    State('protein-fig-holder', 'children'),
    State('protein-datatable-holder','children'),
    State('protein-fig', 'figure'),
    State('protein-list-df-holder', 'children')
)(create_protein_figure_and_table)

app.callback(
    Output('cutoff-value-holder', 'children'),
    Input('tot-intensity-cutoff', 'value'),
    Input('tot-spc-cutoff', 'value'),
    Input('nbr-of-peptides-cutoff', 'value'),
    Input('peptide-intensity-cutoff', 'value'),
    Input('peptide-spc-cutoff', 'value'),
    Input('RT-CCS-checkbox', 'value'),
    Input('proteins-present-in-all-samples-checkbox', 'value')
)(set_cutoffs)

app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])(display_page)

app.callback(
    Output("modal-file", "is_open"),
    [Input("open-modal-file", "n_clicks"), Input("close-modal-file", "n_clicks")],
    [State("modal-file", "is_open")],
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
    Output("modal-normalization", "is_open"),
    [Input("open-modal-normalization", "n_clicks"), Input("close-modal-normalization", "n_clicks")],
    [State("modal-normalization", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-FAQ", "is_open"),
    [Input("open-modal-FAQ", "n_clicks"), Input("close-modal-FAQ", "n_clicks")],
    [State("modal-FAQ", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-feedback", "is_open"),
    [Input("open-modal-feedback", "n_clicks"), Input("close-modal-feedback", "n_clicks")],
    [State("modal-feedback", "is_open")],
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
