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
from methods import amino_acid_piecharts, common_family, all_sample_bar_chart, create_peptide_list_from_trivname
from methods import apply_protein_cutoffs, apply_peptide_cutoffs, get_unique_and_common_proteins
from methods import proteins_present_in_all_samples, create_protein_datatable, create_peptide_datatable, log_intensity, normalize_data, create_length_histogram



app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(id='page-content', fluid=True, className='vh-100'),
 ])

#---------------------------------------PAGE-ELEMENTS------------------------------------------------
file_columns = ['Sample', 'File']


modal_file = html.Div([
    dbc.Button("Files", id="open-modal-file", color='info', className="mr-1"),
    dbc.Tooltip(
        "Import files",
        target='open-modal-file',
        placement='bottom'
    ),
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
                    dbc.Button("Upload files", color = 'primary', id="close-modal-file", className="ml-auto", n_clicks_timestamp=0)
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
                dbc.ModalBody([
                    dbc.Row(dbc.Button('Download protein table', id='download-protein-info'), justify='center', className='mr-1'),
                    dbc.Row(dbc.Button('Download peptide table', id='download-peptide-table'), justify='center', className='mr-1'),
                ]),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-export-data", className="ml-auto")
                ),
            ],
            id="modal-export-data",
            centered=True,
              )
])

protein_tab = dbc.Form([
                    dbc.FormGroup([
                    dbc.Label("Total intensity", className="mr-2"),
                    dbc.Input(placeholder='0', type='number', className='ml-auto', min=0, id='tot-intensity-cutoff', value=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Total spectral count', className='mr-2'),
                    dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, id='tot-spc-cutoff', value=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Number of peptides', className='mr-2'),
                    dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, id='nbr-of-peptides-cutoff', value=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Checklist(
                         options=[
                            {"label": "Only show proteins present in all samples", "value": 'present-in-all-samples'},
                            ],
                            value=[],
                            id='proteins-present-in-all-samples-checkbox')],
                     className="mr-3",)
                ],
                inline=True),

peptide_tab = dbc.Form([
                    dbc.FormGroup([
                    dbc.Label("Intensity", className="mr-2"),
                    dbc.Input(placeholder='0', type='number', className='ml-auto', min=0, id='peptide-intensity-cutoff', value=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Label('Spectral count', className='mr-2'),
                    dbc.Input(placeholder='0', type="number", className='ml-auto', min=0, id='peptide-spc-cutoff', value=0),
                    ], className="mr-3",),
                    dbc.FormGroup([
                    dbc.Checklist(
                         options=[
                            {"label": "RT", "value": 'RT'},
                            {"label": "CSS", "value": 'CSS'},
                            ],
                            value=[],
                            id='RT-CSS-checkbox')],
                     className="mr-3",)
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
            size='m',
            centered=True,
              )

modal_FAQ = html.Div([
    dbc.Button("FAQ", id="open-modal-FAQ", color='info', className="mr-1"),
    dbc.Modal([
                dbc.ModalHeader("FAQ", className="font-weight-bold"),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-FAQ", className="ml-auto")
                ),
            ],
            id="modal-FAQ",
            size='m',
            centered=True,
              )])

modal_feedback = html.Div([
    dbc.Button("Feedback", id="open-modal-feedback", color='info', className="mr-1"),
    dbc.Modal([
                dbc.ModalHeader("Feedback", className="font-weight-bold"),
                dbc.Row(dbc.Textarea(
                    id='textarea-feedback',
                    placeholder='Send your feedback',
                    style={'width': '80%', 'height': 100},
                    
                        ), 
                        align='center', justify='center'),
                
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-modal-feedback", className="ml-auto")
                ),
            ],
            id="modal-feedback",
            size='m',
            centered=True,
              )])

modal_normalization = dbc.Modal([
                dbc.ModalHeader("Normalize data", className="font-weight-bold"),
                dbc.ModalBody([
                    dbc.FormGroup([
                        dbc.Label('log10 area', className='ml-auto'),
                        dbc.Checkbox(id='log-checkbox', checked=True)
                    ]),
                    dbc.FormGroup([
                        dbc.Label('Normalize data', className='ml-auto'),
                        dbc.RadioItems(
                            options=[
                            {'label': 'Normalize on global intensity', 'value': 'global-intensity'},
                            {'label': 'Normalize on housekeeping protein', 'value': 'housekeeping-protein'}
                            ],
                            value='',
                            id='normalization-radioitems',
                        ),
                        dbc.Input(
                            id='housekeeping-protein-input',
                            placeholder='Search protein...',
                            name = 'text',
                            debounce=True,
                            inputMode='latin',
                            minLength=0, maxLength=30,
                            size = '10',
                            list = 'protein-list',
                            className="ml-auto",
                            disabled=True
                        )
                    ]),
                ]),
                dbc.ModalFooter(
                    dbc.Button("Apply", id="close-modal-normalization", className="ml-auto")
                ),
            ],
            id="modal-normalization",
            size='m',
            centered=True,
              )

navbar = dbc.Navbar(
    [
        dbc.NavbarBrand("Eriks och Simons kandidatarbete"),
        modal_file,
        modal_export_data,
        dbc.DropdownMenu(label="Settings",
            children=[
                dbc.DropdownMenuItem("Cutoffs", id="open-modal-cutoff"),
                dbc.DropdownMenuItem("Normalization", id="open-modal-normalization"),
                modal_cutoff,
                modal_normalization,
            ]
        ),
        dbc.Nav([
        modal_FAQ,
        modal_feedback,
        ],
        navbar=True,
        className="ml-auto",)
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
            dbc.CardHeader("How to use!"),
            dbc.CardBody("This will display steps and links on how to use the app."),
            ]),
            id="how-to-use-collapse",
            className="mb-4"
        ),
        dbc.Tooltip(
            "View guideline on how to use (NAMN)",
            target="how-to-use-collapse-button", 
            placement="right"
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
        dbc.Tooltip(
            "View samples",
            target="sample-collapse-button",
            placement="right",
            #style={},
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

protein_fig_radioitems = html.Div([
    dbc.Label("Select difference metric"),
    dbc.RadioItems(
        options=[
        {'label': 'Area', 'value': 'area'},
        {'label': 'Spectral Count', 'value': 'spectral_count'}
        ],
        value='area',
        id='protein-radioitems',
        inline=True,
    )
])

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
    dbc.Label("Select difference metric"),
    dbc.RadioItems(
        options=[
        {'label': 'Area', 'value': 'area'},
        {'label': 'Spectral Count', 'value': 'spectral_count'}
        ],
        value='area',
        id='peptide-radioitems',
        inline=True,
    )
])

peptide_fig_radioitems_sum_or_mean = html.Div([
    dbc.Label("Sum or mean"),
    dbc.RadioItems(
        options=[
        {'label': 'sum', 'value': False},
        {'label': 'mean', 'value': True}
        ],
        value=False,
        id='sum-or-mean-radio',
        inline=True,
    )
])

peptide_fig = html.Div([
        html.H3('Peptide View'),
        dbc.Row([
            dbc.Col(peptide_fig_radioitems_sum_or_mean),
            dbc.Col(peptide_fig_radioitems)    ,
            dbc.Col(dbc.Button('Choose a protein', disabled=True, id='generate-peptide-fig', color='success'))             
        ]),
        dcc.Loading(type='cube', color = '#76b382',
            children=dcc.Graph(id='peptide-fig', figure={}, config={'displaylogo': False})
        )
        ])

amino_acid_radioitems = html.Div([
        dbc.Label("Select difference metric"), 
        dbc.RadioItems(
            options=[
                {"label": "Area", "value": 'area'},
                {"label": "Spectral Count", "value": 'spectral_count'}, 
            ], 
            value='area',
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
                    dcc.Graph(id='complete-aa-seq-fig-g1', figure={}, config={'displaylogo': False})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='first-aa-fig-g1', figure={}, config={'displaylogo': False})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='last-aa-fig-g1', figure={}, config={'displaylogo': False})], width={'size':3}),
            ]),
            dbc.Row([
                dbc.Col([
                    dcc.Graph(id='complete-aa-seq-fig-g2', figure={}, config={'displaylogo': False})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='first-aa-fig-g2', figure={}, config={'displaylogo': False})], width={'size':3}),
                dbc.Col([
                    dcc.Graph(id='last-aa-fig-g2', figure={}, config={'displaylogo': False})], width={'size':3}),
            ])]
        )
    ])

peptide_length_dropdown = dcc.Dropdown(
    id= 'peptide-length-dropdown',
    placeholder='Select view',
    value='',
    options=[
        {'label': 'Complete proteome', 'value': 'complete-proteome-length'},
        {'label': 'Selected protein', 'value': 'selected-protein-length'},
    ],
    
)    

peptide_length_figs = html.Div([
    html.H3('Peptide Length'),
    dbc.Row([
        dbc.Col(peptide_length_dropdown, width={'size':2})
    ]),
    dcc.Loading(type='cube', color = '#76b382',
        children=[    
        dbc.Row([
            dbc.Col([
                dcc.Graph(id='peptide-length-fig-g1', figure={}, config={'displaylogo': False})], width={'size':5}),
            dbc.Col([
                dcc.Graph(id='peptide-length-fig-g2', figure={}, config={'displaylogo':False})], width={'size':5}),
        ])]
    )
])

protein_info = html.Div(dash_table.DataTable(
            id='protein-info-table',
            sort_action='native',
            fixed_rows={'headers': True},
            #filter_action='native',
            virtualization=True,
            row_selectable="multi",
            export_format='xlsx',
            selected_rows=[],
            css=[{'selector':'.export','rule':'position:font-type:Roboto;color:black;background-color:#FAFAFA;border-color:#FAFAFA;border:1px solid transparent'}],
            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : 'rgb(182, 224, 194)'
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
            sort_action='native',
            fixed_rows={'headers': True},
            #filter_action='native',
            virtualization=True,
            row_selectable="multi",
            export_format='xlsx',
            selected_rows=[],
            css=[{'selector':'.export','rule':'position:font-type:Roboto;color:black;background-color:#FAFAFA;border-color:#FAFAFA;border:1px solid transparent'}],
            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : 'rgb(182, 224, 194)'
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


hidden_divs = html.Div([
    html.Div(id='cutoff-value-holder', style={'display': 'none'}),
    html.Div(id='protein-list-df-holder', style={'display': 'none'}),
    html.Div(id='peptide-list-df-holder', style={'display': 'none'}),
    html.Div(id='df_g1-holder', style={'display':'none'}),
    html.Div(id='df_g2-holder', style={'display':'none'}),
    html.Div(id='protein-datatable-holder', style={'display':'none'}),
    html.Div(id='peptide-data-holder', style={'display':'none'}),
    html.Div(id='protein-fig-holder', style = {'display':'none'}),
    html.Div(id='normalization-holder', style = {'display':'none'}),
    html.Div(id='housekeeping-protein-holder', style = {'display':'none'})
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
    ]),

    dbc.Row([
        dbc.Col(peptide_length_figs),
    ]),  
    hidden_divs,
], fluid=True)


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


def update_file_list(contents, filename):
    if filename:
        file_list = []
        i=0
        for f in filename:
            s = [f"S{i}", f]
            file_list.append(s)
            i+=1  
        master_df = update_data_frame(contents, filename)
        df = pd.DataFrame(file_list, columns = ['Sample', 'File'])
        return df.to_dict('rows'), df.to_dict('rows'), master_df.to_json()
    else:
        return [],[], []

def update_data_frame(contents, filename):
    decoded_list = []
    for f, content in zip(filename, contents):
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        decoded_list.append(io.BytesIO(decoded))
    dfs = make_peptide_dfs(decoded_list)
    master_df = concatenate_dataframes(dfs)
    master_df = log_intensity(master_df)
    return master_df


def set_cutoffs(tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT_CSS_checkbox, proteins_present_in_all_samples):
    RT = False
    CSS = False
    present_in_all_samples = False
    if 'RT' in RT_CSS_checkbox:
        RT=True
    if 'CSS' in RT_CSS_checkbox:
        CSS=True
    if 'present-in-all-samples' in proteins_present_in_all_samples:
        proteins_present_in_all_samples = True
    return [tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples]


def process_protein_data(apply_normalization_n_clicks, n_clicks_close_file, apply_cutoffs_button, cutoff_values, df_g1, df_g2, radioitems_normalization, housekeeping_protein):
    if apply_cutoffs_button:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = cutoff_values
    else:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = 0,0,0,0,0,False,False, False
    triv_names = []
    protein_fig = {}
    protein_list = []
    protein_list_cutoff = []
    if n_clicks_close_file and df_g1 and df_g2:
        g1 = pd.read_json(df_g1)
        g2 = pd.read_json(df_g2)
        master = merge_dataframes(g1,g2)
        protein_list = create_protein_list(master)
        protein_list_cutoff = apply_peptide_cutoffs(protein_list, area=pep_intensity_co, spc=pep_spc_co, rt=RT, css=CSS)
        protein_list_cutoff = apply_protein_cutoffs(protein_list_cutoff, nbr_of_peptides=nbr_of_peptides_co, tot_area=tot_intensity_co, tot_spc=tot_spc_co)
        
        if present_in_all_samples:
            protein_list_cutoff = proteins_present_in_all_samples(protein_list_cutoff)
        if apply_normalization_n_clicks:
            if 'global-intensity' in radioitems_normalization:
                protein_list_cutoff = normalize_data(protein_list_cutoff, housekeeping_protein=False)
            elif 'housekeeping-protein' in radioitems_normalization and housekeeping_protein != '':
                protein_list_cutoff = normalize_data(protein_list_cutoff, housekeeping_protein = housekeeping_protein)
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
    print(changed_id)
    highlighted_triv_names = 'Choose protein'
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
        if 'area' in protein_radioitems_value:
            difference_metric = 'area'
            columns = ['Protein','UniProt id','#peptides g1','#peptides g2', 'intensity_g1','intensity_g2', 'Protein family','p-value_area']
            sort = ['intensity_g1', 'intensity_g2']
            x_label = 'Group 1 log(intensity)'
            y_label = 'Group 2 log(intensity)'
        else:
            difference_metric = 'spc'
            columns = ['Protein','UniProt id','#peptides g1','#peptides g2','spc_g1','spc_g2', 'Protein family','p-value_spc']
            sort = ['spc_g1','spc_g2']
            x_label = 'Group 1 spectral count'
            y_label = 'Group 2 spectral count'
        df_protein_info.sort_values(by=sort, ascending=False, inplace=True)
        protein_info_data = df_protein_info.to_dict('rows')
        protein_info_columns=[{"name": str(i), "id": str(i)} for i in columns]
        if str(changed_id) == 'generate-protein-graph.n_clicks' or str(changed_id) == 'protein-radioitems.value' or str(changed_id) == 'protein-checkbox.value':
            if checkbox_values and 'show-stdev' in checkbox_values and 'show-pfam' in checkbox_values:
                protein_fig = create_protein_fig(df_fig, protein_list, show_pfam=True, show_stdev = True,  difference_metric = difference_metric)
            elif checkbox_values and 'show-stdev' in checkbox_values:
                protein_fig = create_protein_fig(df_fig, protein_list, show_stdev = True, difference_metric = difference_metric)
            elif checkbox_values and 'show-pfam' in checkbox_values:
                protein_fig = create_protein_fig(df_fig, protein_list, show_pfam=True, difference_metric = difference_metric)
            else:
                protein_fig = create_protein_fig(df_fig, protein_list, difference_metric = difference_metric)
        return protein_fig, protein_info_data, protein_info_columns, disabled, str(highlighted_triv_names[0])
    else:
        return {}, [], [], disabled, 'Choose protein'

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
        columns = ['Peptide','Start','End','Intensity_g1','Intensity_g2']
        sort = ['Intensity_g1', 'Intensity_g2']
    else:
        columns = ['Peptide','Start','End','spc_g1','spc_g2']
        sort = ['spc_g1', 'spc_g2']

    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if changed_id == 'generate-peptide-fig.n_clicks':
        protein_list = json_to_protein_list(protein_list_json)
        trivname = button_label.split(' ')[-1]
        peptide_list = create_peptide_list_from_trivname(protein_list, str(trivname))
        peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=sum_or_mean_radio, difference_metric=peptide_radioitems_value)
        df_peptide_info = create_peptide_datatable(peptide_list)
        df_peptide_info.fillna(0, inplace=True)
        df_peptide_info.sort_values(by=sort, ascending=False, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows')
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in columns]
        return peptide_fig,  peptide_table_data, peptide_table_columns, peptide_list_to_json(peptide_list)
    elif peptide_list_json and (changed_id == 'sum-or-mean-radio.value' or 'peptide-radioitems.value'):
        peptide_list = json_to_peptide_list(peptide_list_json)
        peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=sum_or_mean_radio, difference_metric=peptide_radioitems_value)
        df_peptide_info = create_peptide_datatable(peptide_list)
        df_peptide_info.fillna(0, inplace=True)
        df_peptide_info.sort_values(by=sort, ascending=False, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows')
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in columns]
        return peptide_fig, peptide_table_data, peptide_table_columns, peptide_list_json
    
    else:
        return {}, [], [], []


def amino_acid_dropdown(n_clicks_complete_proteome, n_clicks_selected_protein, radioitem_value, protein_list_json, peptide_list_json):
    
    if n_clicks_complete_proteome > n_clicks_selected_protein and protein_list_json:
        protein_list = json_to_protein_list(protein_list_json)
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(protein_list, peptide_or_protein_list = 'protein_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Complete proteome'
    elif n_clicks_selected_protein > n_clicks_complete_proteome and peptide_list_json:
        peptide_list = json_to_peptide_list(peptide_list_json)
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_list, peptide_or_protein_list = 'peptide_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Selected Protein'
    else:
        return  {},{},{}, {},{},{}, ''


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

def peptide_length_dropdown(length_dropdown_values, protein_list_json, peptide_list_json):
    
    if 'complete-proteome-length'  in length_dropdown_values and protein_list_json:
        protein_list = json_to_protein_list(protein_list_json)
        length_fig_g1, length_fig_g2 = create_length_histogram(protein_list, peptide_or_protein_list='protein_list')
        return length_fig_g1, length_fig_g2
    elif 'selected-protein-length' in length_dropdown_values and peptide_list_json:
        peptide_list = json_to_peptide_list(peptide_list_json)
        length_fig_g1, length_fig_g2 = create_length_histogram(peptide_list, peptide_or_protein_list='peptide_list')
        return length_fig_g1, length_fig_g2
    else:
        return {}, {}


app.callback(
    Output('generate-protein-graph', 'disabled'),
    Input('protein-list', 'children')
)(enable_input_search_protein)

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
    Output('complete-aa-seq-fig-g1', 'figure'),
    Output('first-aa-fig-g1', 'figure'),
    Output('last-aa-fig-g1', 'figure'),
    Output('complete-aa-seq-fig-g2', 'figure'),
    Output('first-aa-fig-g2', 'figure'),
    Output('last-aa-fig-g2', 'figure'),
    Output('complete-or-selected', 'children'),
    Input('complete-proteome', 'n_clicks_timestamp'),
    Input('selected-protein','n_clicks_timestamp'),
    Input('aa-radioitems', 'value'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(amino_acid_dropdown)

app.callback(
    Output('peptide-length-fig-g1', 'figure'),
    Output('peptide-length-fig-g2', 'figure'),
    Input('peptide-length-dropdown', 'value'),
    State('protein-list-df-holder', 'children'),
    State('peptide-list-df-holder', 'children'),
)(peptide_length_dropdown)


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
    [Input("upload-data-1", "contents"), Input("upload-data-1", "filename")],
)(update_file_list)

app.callback(
    Output("output-filename-2", "data"),
    Output("sample-collapse-2", "data"),
    Output('df_g2-holder','children'),
    [Input("upload-data-2", "contents"), Input("upload-data-2", "filename")],
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
    Input('RT-CSS-checkbox', 'value'),
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
