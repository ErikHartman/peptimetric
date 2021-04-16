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

from methods import make_peptide_dfs, concatenate_dataframes, merge_dataframes, create_protein_list, protein_graphic_plotly, create_peptide_list, stacked_samples_peptide
from methods import amino_acid_piecharts, common_family, all_sample_bar_chart
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
            dbc.Col(protein_fig_radioitems)
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

peptide_fig = html.Div([
        html.H3('Peptide View'),
        dbc.Row([
            dbc.Col(
            dbc.DropdownMenu(label='peptide dropdown',
            children = [
                dbc.DropdownMenuItem("Sum", id="peptide-dropdown-sum", n_clicks_timestamp=0),
                dbc.DropdownMenuItem("Mean", id="peptide-dropdown-mean", n_clicks_timestamp=0),
            ])),
            dbc.Col(
                peptide_fig_radioitems
            )                 
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


peptide_info = html.Div([dash_table.DataTable(
            id='peptide-info-table',
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
)]),



hidden_divs = html.Div([
    html.Div(id='cutoff-value-holder', style={'display': 'none'}),
    html.Div(id='protein-lists-holder', style={'display': 'none'}),
    html.Div(id='peptide-lists-holder', style={'display': 'none'}),
    html.Div(id='df_g1-holder', style={'display':'none'}),
    html.Div(id='df_g2-holder', style={'display':'none'}),
    html.Div(id='protein-data-holder', style={'display':'none'}),
    html.Div(id='peptide-data-holder', style={'display':'none'}),
    html.Div(id='trivial-name-holder', style={'display':'none'})
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

protein_lists = []

def create_protein_fig(n_clicks, checkbox_values, apply_cutoffs_button, protein_radioitems_value, normalization_radioitems, rows, selected_rows, housekeeping_protein, cutoff_values, df_g1, df_g2, protein_fig):
    
    if selected_rows is not None:
        df= pd.DataFrame(selected_rows)
        print(df)

    if apply_cutoffs_button:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = cutoff_values

    else:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = 0,0,0,0,0,False,False, False
    triv_names = []
    protein_fig = {}
    protein_list = []
    protein_list_cutoff = []
    if n_clicks and df_g1 and df_g2:
        g1 = pd.read_json(df_g1)
        g2 = pd.read_json(df_g2)
        master = merge_dataframes(g1,g2)
        protein_list = create_protein_list(master)
        protein_lists.append(protein_list)
        protein_list_cutoff = apply_peptide_cutoffs(protein_list, area=pep_intensity_co, spc=pep_spc_co, rt=RT, css=CSS)
        protein_list_cutoff = apply_protein_cutoffs(protein_list_cutoff, nbr_of_peptides=nbr_of_peptides_co, tot_area=tot_intensity_co, tot_spc=tot_spc_co)
        if present_in_all_samples:
            protein_list_cutoff = proteins_present_in_all_samples(protein_list_cutoff)
        if 'global-intensity' in normalization_radioitems:
            protein_list_cutoff = normalize_data(protein_list_cutoff, housekeeping_protein =False)
        if 'housekeeping-protein' in normalization_radioitems and housekeeping_protein != '':
            protein_list_cutoff = normalize_data(protein_list_cutoff, housekeeping_protein = housekeeping_protein)
        unique_protein_list, common_protein_list = get_unique_and_common_proteins(protein_list_cutoff)
        df_protein_info = create_protein_datatable(protein_list_cutoff, protein_radioitems_value)
        df_protein_info.fillna(0, inplace=True)
        protein_info_data = df_protein_info.to_dict('rows')
        protein_info_columns=[{"name": str(i), "id": str(i)} for i in df_protein_info.columns]
            
        if len(common_protein_list) > 1:
            for protein in protein_list_cutoff:
                triv_names.append(html.Option(value=protein.get_trivial_name()))
    if df_g1 and df_g2:
        if checkbox_values and 'show-stdev' in checkbox_values and 'show-pfam' in checkbox_values:
            protein_fig = protein_graphic_plotly(common_protein_list, difference_metric= protein_radioitems_value, show_pfam=True, show_stdev = True)
        elif checkbox_values and 'show-stdev' in checkbox_values:
            protein_fig = protein_graphic_plotly(common_protein_list, difference_metric=protein_radioitems_value, show_stdev = True)
        elif checkbox_values and 'show-pfam' in checkbox_values:
            protein_fig = protein_graphic_plotly(common_protein_list, difference_metric=protein_radioitems_value, show_pfam = True)
        else:

            protein_fig = protein_graphic_plotly(common_protein_list, difference_metric=protein_radioitems_value)
        return protein_fig, triv_names, protein_info_data, protein_info_columns, df_protein_info.to_json()

    
    else:
        return {}, [], [], [], []

peptide_lists=[]
def create_peptide_fig(clickData, search_protein, n_clicks_sum, n_clicks_mean, cutoff_values, apply_cutoffs_button, peptide_radioitems_value,):
    protein_accession = ''
    search_text = ''
    if apply_cutoffs_button:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = cutoff_values
    else:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CSS, present_in_all_samples = 0,0,0,0,0,False,False, False

    if len(protein_lists) > 0:
        protein_list_cutoff = apply_peptide_cutoffs(protein_lists[-1], area=pep_intensity_co, spc=pep_spc_co, rt=RT, css=CSS)
        protein_list_cutoff = apply_protein_cutoffs(protein_list_cutoff, nbr_of_peptides=nbr_of_peptides_co, tot_area=tot_intensity_co, tot_spc=tot_spc_co)
    if search_protein != '' and len(protein_lists) > 0:
        for protein in protein_list_cutoff:
            if search_protein == protein.get_trivial_name():
                protein_accession = protein.get_id()
                search_text = ''
    elif clickData:
        protein_accession = clickData['points'][0]['customdata'][-1]
        search_text=''

    if protein_accession != '':
        peptide_list = create_peptide_list(protein_list_cutoff, str(protein_accession))
        peptide_lists.append(peptide_list)
        if n_clicks_sum == 0 and n_clicks_mean == 0 or n_clicks_sum > n_clicks_mean:
            peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=False, difference_metric=peptide_radioitems_value)
        else:
            peptide_fig = stacked_samples_peptide(peptide_list, show_difference='show', show_weight ='show', average=True, difference_metric=peptide_radioitems_value)
        df_peptide_info = create_peptide_datatable(peptide_list, peptide_radioitems_value)
        df_peptide_info.fillna(0, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows'),
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in df_peptide_info.columns],
        return peptide_fig, search_text,  peptide_table_data, peptide_table_columns, df_peptide_info.to_json()
    
    else:
        return {}, search_text, [], [], []


def amino_acid_dropdown(n_clicks_complete_proteome, n_clicks_selected_protein, radioitem_value):
    if n_clicks_complete_proteome > n_clicks_selected_protein and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(protein_lists[-1], peptide_or_protein_list = 'protein_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Complete proteome'
    elif n_clicks_selected_protein > n_clicks_complete_proteome and len(protein_lists)>0:
        complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_lists[-1], peptide_or_protein_list = 'peptide_list', difference_metric = radioitem_value)
        return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2, 'Selected Protein'
    else:
        return  {},{},{}, {},{},{}, ''


def generate_hover_graphs(hoverData, protein_radioitems_value):
    if hoverData:
        accession = hoverData['points'][0]['customdata'][-1]
        protein_list = protein_lists[-1]
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

def highlight_protein_fig(highlighted_triv_names, triv_names_holder, protein_fig):
    if highlighted_triv_names:
        highlighted_triv_names = clickData['points'][0]['customdata'][-1]
    marker_color_list = ['rgba(0,0,0,0)' for n in range(len(triv_names_holder))]
    for i in range(len(triv_names_holder)):
        if highlighted_triv_names == str(triv_names_holder[i]):
            marker_color_list[i] = red['medium']
            protein_fig.update_traces(marker=dict(line=dict(width=2, color=marker_color_list)),
                selector=dict(mode='markers'))
            return protein_fig



def peptide_length_dropdown(length_dropdown_values):
    if 'complete-proteome-length'  in length_dropdown_values and len(protein_lists) > 0:
        length_fig_g1, length_fig_g2 = create_length_histogram(protein_lists[-1], peptide_or_protein_list='protein_list')
        return length_fig_g1, length_fig_g2
    elif 'selected-protein-length' in length_dropdown_values and len(protein_lists) > 0:
        length_fig_g1, length_fig_g2 = create_length_histogram(peptide_lists[-1], peptide_or_protein_list='peptide_list')
        return length_fig_g1, length_fig_g2
    else:
        return {}, {}


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
    Output('peptide-length-fig-g1', 'figure'),
    Output('peptide-length-fig-g2', 'figure'),
    Input('peptide-length-dropdown', 'value'),
)(peptide_length_dropdown)

app.callback(
    Output('peptide-fig', 'figure'),
    Output('search-protein', 'value'),
    Output('peptide-info-table', 'data'),
    Output('peptide-info-table', 'columns'),
    Output('peptide-data-holder', 'children'),
    Input('protein-fig', 'clickData'),
    Input('search-protein', 'value'),
    Input('peptide-dropdown-sum', 'n_clicks_timestamp'),
    Input('peptide-dropdown-mean','n_clicks_timestamp'),
    Input('cutoff-value-holder','children'),
    Input('close-modal-cutoff', 'n_clicks'),
    Input('peptide-radioitems', 'value')
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
    Output('protein-fig', 'figure'),
    Output('protein-list', 'children'),
    Output('protein-info-table','data'),
    Output('protein-info-table','columns'),
    Output('protein-data-holder', 'children'),
    Input("close-modal-file", "n_clicks_timestamp"),
    Input('protein-checklist','value'),
    Input('close-modal-cutoff', 'n_clicks'),
    Input('protein-radioitems', 'value'),
    Input('normalization-radioitems', 'value'),
    Input('housekeeping-protein-input','value'),
    Input('protein-info-table', "derived_virtual_data"),
    Input('protein-info-table', "derived_virtual_selected_rows"),
    State('cutoff-value-holder', 'children'),
    State('df_g1-holder', 'children'),
    State('df_g2-holder', 'children'),
    State('protein-fig','figure'),
    )(create_protein_fig)

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
