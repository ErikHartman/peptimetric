
import base64
import io
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go
from dash import no_update

from os import listdir

from methods import pre_process_peptide_fig
from methods import make_peptide_dfs, concatenate_dataframes, merge_dataframes, create_protein_df_fig, create_protein_fig , create_peptide_fig
from methods import amino_acid_piecharts, all_sample_bar_chart, protein_create_protein_list
from methods import apply_protein_cutoffs, apply_peptide_cutoffs, create_venn_bar
from methods import proteins_present_in_all_samples, create_protein_datatable, create_peptide_datatable, log_intensity, normalize_data, create_length_histogram
from texts_for_webapp import how_to_use, Documentation, contact_text
from dash_extensions.enrich import Dash, ServersideOutput, Output, Input, State


app = Dash(__name__, title='peptimetric', external_stylesheets=[dbc.themes.SANDSTONE], suppress_callback_exceptions=True)
server=app.server

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(id='page-content', fluid=True, className='vh-100'),
 ])

#---------------------------------------PAGE-ELEMENTS------------------------------------------------
file_columns = ['Sample', 'File']
sample_files = pd.read_csv('./example-files/all-files.csv')
sample_files = protein_create_protein_list(sample_files, 'homo-sapiens')

modal_file = html.Div([

    dbc.Button("Upload files", id="open-modal-file", color='secondary',  outline=True, className="mr-1"),
        dbc.Modal([
                dbc.ModalHeader("Upload files", className="font-weight-bold"),
                dbc.ModalBody([
                    dbc.FormGroup([
                            dbc.Label("Select species", html_for="select-species", width=3, style={'padding-left':20, 'padding-right':5}),
                            dbc.Col(dcc.Dropdown(
                                        id= 'select-species',
                                        placeholder='Select species...',
                                        value='homo-sapiens',
                                        options=[
                                            {'label': 'Homo sapiens', 'value': 'homo-sapiens'},
                                            {'label': 'Sus scrofa (Pig)', 'value': 'pig'},
                                            {'label': 'Rattus norvegicus (Rat)', 'value': 'rat'},
                                            {'label': 'Cricetulus griseus (Hamster)', 'value': 'hamster'},
                                            {'label': 'Mus musculus (Mouse)', 'value': 'mouse'},
                                            {'label': 'Danio rerio (Zebra fish)', 'value': 'zebra-fish'},
                                            {'label': 'Drosophila melanogaster', 'value': 'drosophila'},
                                            {'label': 'Caenorhabditis elegans', 'value': 'c-elegans'},
                                            {'label': 'Candida albicans', 'value': 'candida'},
                                            {'label': 'Escherichia coli', 'value': 'ecoli'},

                                        ],
                                        
                                    )   ,width=7,style={'padding-left':6}),
                        ],
                        row=True,
                               ),
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
                                        'backgroundColor' : '#a4d694'
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
                                        'backgroundColor' : '#a4d694'
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
                            
                    ],),
                    dbc.Label('*upload 3 or more files for statistical analysis'),
                    ]),
                dbc.ModalFooter([
                    dbc.Button("Close", color = 'secondary', id="close-modal-file-2", outline=True, className="mr-auto", n_clicks_timestamp=0),
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
                dbc.Tabs(className='custom-tabs',
                    children=[
                    dbc.Tab(peptide_tab, label='Peptide', ),
                    dbc.Tab(protein_tab, label='Protein'),
                    
                ]),
                dbc.Label('*cutoffs will be applied after normalization', style={'padding':10}),
                dbc.ModalFooter([
                    dbc.Button("Close", id="close-modal-cutoff-2", color="secondary", outline=True, className="mr-auto"),
                    dbc.Button("Apply", id="close-modal-cutoff", color="primary", className="ml-auto")
                ]),
            ],
            id="modal-cutoff",
            size='m',
            centered=True,
              )



modal_feedback = html.Div([
    dbc.Button("Feedback", id="open-modal-feedback", color='secondary', outline=True, className="mr-1"),
    dbc.Modal([
        dbc.ModalHeader('Feedback'),
                dbc.ModalBody(
                    html.Div(
                        html.P(contact_text)
                    )
                ),
                
                dbc.ModalFooter([
                    dbc.Button("Close", id="close-modal-feedback", color='secondary', outline=True,className="mr-auto")
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
                                    {'label': 'Normalize on global values', 'value': 'global-intensity'},
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
                
                dbc.ModalFooter([
                    dbc.Button("Close", id="close-modal-normalization-2", color="secondary", outline=True, className="mr-auto"),
                    dbc.Button("Apply", id="close-modal-normalization", color="primary", className="ml-auto")
                ]),
            ],
            id="modal-normalization",
            size='m',
            centered=True,
              )


navbar = dbc.Navbar(
    [
        dbc.NavLink("Peptimetric", href = '/', style = {'color':'grey', 'font-size':20,  'font-weight':'bold', 'font':'Roboto', 'text-transform':'lowercase'} ),
        modal_file,
        dbc.Button('Normalization', id="open-modal-normalization", color='secondary', outline=True, className='mr-1'),
        normalization_modal,
        dbc.Button("Cutoffs", id="open-modal-cutoff", color='secondary', outline=True, className='mr-1'),
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
    dcc.Loading(color = '#cf597e', style={'backgroundColor': 'transparent'}, className = 'loader-wrapper', fullscreen = True, type = 'default', id='process-data-loading', children = [html.Datalist( id = 'protein-list', children=[])]),
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
            dbc.CardFooter(dbc.Button('Load example data', id='load-sample-data', color='primary', n_clicks_timestamp=0)
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
                                        'backgroundColor' : '#a4d694'
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
                                        'backgroundColor' : '#a4d694'
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
        placeholder = 'Select abundance metric...',
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
        html.Img(src = app.get_asset_url ('scatter.jpg'), style={'height':'3%', 'width':'3%'}),
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
        dcc.Loading(type='cube', color = '#a4d694',
            children=dcc.Graph(id='protein-fig', figure={}, config={'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
        },'displaylogo':False})
        )
        
        ])

all_samples_protein_fig = html.Div([
    dcc.Graph(id='hover-all-protein-samples', figure={}, style={'height': 300, 'width':500},
    config={
        'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
        },
        'displaylogo': False, 'modeBarButtonsToRemove': ['toggleSpikelines','hoverCompareCartesian','zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d','hoverClosestGl2d',
                'hoverClosestGl2d','hoverClosestPie', 'hoverClosestCartesian', 'autoScale2d', 'resetScale2d']})])

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
        html.Img(src = app.get_asset_url ('bar.jpg'), style={'height':'3%', 'width':'3%'}),
        html.H3('Peptide View'),
    ]),
        dbc.Row([
            dbc.Col(peptide_fig_radioitems_sum_or_mean),
            dbc.Col(peptide_fig_radioitems)    ,
            dbc.Col(dbc.Button('Choose a protein', disabled=True, id='generate-peptide-fig', color='success'))             
        ]),
        dcc.Loading(type='cube', color = '#a4d694',
            children=dcc.Graph(id='peptide-fig', figure={}, config={'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2
        },'displaylogo': False})
        ),
        ])

amino_acid_radioitems = html.Div([
        dbc.Label('Abundance metric'),
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
        dcc.Loading(type='cube', color = '#a4d694',
            children=[ dbc.Row([
                dbc.Col(dcc.Graph(id='aa-fig', figure={}, style={'height': '700px'}, config={'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2 
        },'displaylogo': False, 'modeBarButtonsToRemove': ['toggleSpikelines','hoverCompareCartesian','zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d','hoverClosestGl2d',
                'hoverClosestPie', 'hoverClosestCartesian', 'autoScale2d', 'resetScale2d']})),
            ]) 
            ]
        )
    ])

peptide_length_fig = html.Div([
    html.H3('Peptide Length'),
    dcc.Loading(type='cube', color = '#a4d694',
        children=[    
        dbc.Row([
            dbc.Col(dcc.Graph(id='peptide-length-fig', figure={}, config={'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2
        },'displaylogo': False, 'modeBarButtonsToRemove': ['toggleSpikelines','hoverCompareCartesian','zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d','hoverClosestGl2d',
                'hoverClosestGl2d','hoverClosestPie', 'hoverClosestCartesian', 'autoScale2d', 'resetScale2d']}))
        ]),])
])

venn_bar_fig = html.Div([
    html.H3('Peptide overlap'),
    dcc.Loading(type='cube', color = '#a4d694',
        children=[    
        dbc.Row([
            dbc.Col([
                dcc.Graph(id='venn-bar', figure={}, config={'toImageButtonOptions': {
            'format': 'svg',
            'scale': 2
        },'displaylogo': False, 'modeBarButtonsToRemove': ['toggleSpikelines','hoverCompareCartesian','zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d','hoverClosestGl2d',
                'hoverClosestGl2d','hoverClosestPie', 'hoverClosestCartesian', 'autoScale2d', 'resetScale2d']})]),
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
            export_format='csv',
            selected_rows=[],
            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : '#a4d694'
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
            row_selectable="multi",
            export_format='csv',
            selected_rows=[],

            style_data_conditional = [{
                'if' : {'row_index':'odd'},
                'backgroundColor' : '#a4d694'
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
    dbc.Row([html.Img(src = app.get_asset_url ('document.jpg'), style={'height':'40px','width':'auto'}),
        (html.H1('Documentation', style={'margin-top':'5px','margin-bottom':'10px'}))]),
    html.Hr(),
    Documentation,
    
])

bottom_navbar = html.Div(dbc.Navbar([
    dbc.Nav([
    dbc.NavLink("Peptimetric", href = '/', style = {'color':'grey', 'font-size':20,  'font-weight':'bold', 'font':'Roboto', 'text-transform':'lowercase'}),
    dbc.NavbarBrand('by: Erik Hartman, Simon Mahdavi', style={'color':'#808080','font-size':14, 'margin-left':15}), 
    ],vertical='md',),
    dbc.Nav([
        dbc.NavLink('Documentation', href='/documentation'),
        dbc.NavLink('Home', href='/')
        ],
        navbar=True,
        horizontal='start',
        vertical='md',
        style={'margin-right':30},
        className="ml-auto")], className='navbar-bottom', fixed=True))


hidden_divs = html.Div([
    dcc.Store(id='cutoff-value-holder'), 
    dcc.Store(id='protein-df'),
    dcc.Store(id='protein-df-cutoff'),
    dcc.Store(id='peptide-df'),
    dcc.Store(id='df_g1-holder'),
    dcc.Store(id='df_g2-holder'),
    dcc.Store(id='peptide-data-holder'),
    dcc.Store(id='protein-df-fig-holder'),
    dcc.Store(id='normalization-holder'),
    dcc.Loading(color = '#cf597e', style={'backgroundColor': 'transparent'}, fullscreen = True, type = 'default', id='process-data-loading', children = [dcc.Store(id='protein-datatable-holder')]),
    dcc.Store(id='housekeeping-protein-holder'),
    dcc.Store(id='processed-peptide-data')
])

hidden_divs_documentation = html.Div([
    hidden_divs,
    sample_collapse,
    search_protein,
    protein_fig_radioitems,
    how_to_use_collapse,
], style={'display':'none'})
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
        html.Img(src = app.get_asset_url ('pie.jpg'), style={'height':'2%', 'width':'2%'}),
        html.H3('General characteristics', style={'bottom-margin':0}),
    ]),
    dbc.Row([dbc.Col(amino_acid_pie_dropdown, width=2, style={'padding':15}),
    dbc.Col(amino_acid_radioitems, style={'padding':15})]),
    dbc.Row([
        dbc.Col(peptide_length_fig, width={'size':8}),
        dbc.Col(venn_bar_fig, width={'size':4}),
    ]),  
    dbc.Row(dbc.Col(amino_acid_figs)),
    
    dbc.Row([dbc.Col(bottom_navbar, width={'size': 12}, style={'padding':0},)]),

    hidden_divs,
], fluid=True, style={'padding':0})

documentation_page = dbc.Container([
     dbc.Row([
        dbc.Col(navbar, width={"size":12}, className="mb-4", style={'padding-right':0, 'padding-left':0})
    ]),
    dbc.Row([
        documentation,
    ]),
    dbc.Row([dbc.Col(bottom_navbar, width={'size': 12}, style={'padding':0},)]),
    hidden_divs_documentation,

    
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

def toggle_modal(n1, n2, n3, is_open):
    if n1 or n2 or n3:
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
        return df.to_dict('rows'), df.to_dict('rows'), master_df
    else:
        return [], [], pd.DataFrame()

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


def apply_cutoffs_to_protein_list(master_df, apply_normalization_n_clicks, apply_cutoffs_button, cutoff_values, radioitems_normalization, housekeeping_protein):
    if cutoff_values:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CCS, present_in_all_samples = cutoff_values
    else:
        tot_intensity_co, tot_spc_co, nbr_of_peptides_co, pep_intensity_co, pep_spc_co, RT, CCS, present_in_all_samples = 0,0,0,0,0,False,False,False
    triv_names = []
    if not master_df.empty:
        if radioitems_normalization and 'global-intensity' in radioitems_normalization:
            master_df = normalize_data(master_df, housekeeping_protein=False)
        elif radioitems_normalization and 'housekeeping-protein' in radioitems_normalization and housekeeping_protein != '':
            master_df = normalize_data(master_df, housekeeping_protein = housekeeping_protein)
        if cutoff_values and cutoff_values != [0,0,0,0,0,False,False,False]:
            master_df = apply_peptide_cutoffs(master_df, area=pep_intensity_co, spc=pep_spc_co, rt=RT, ccs=CCS)
            master_df = apply_protein_cutoffs(master_df, nbr_of_peptides=nbr_of_peptides_co, tot_area=tot_intensity_co, tot_spc=tot_spc_co)
        if len(master_df.index) < 1:
            return [], pd.DataFrame(), []
        if present_in_all_samples:
            master_df = proteins_present_in_all_samples(master_df)
        trivnames = master_df['trivname'].unique()
        for name in trivnames:
            triv_names.append(html.Option(value=name))
        return triv_names, master_df

    else:

        return [], pd.DataFrame()

def make_protein_list(n, load_sample_file_n_clicks, df_g1, df_g2, species):
    if load_sample_file_n_clicks and load_sample_file_n_clicks > n:
        return sample_files
    if n and n > load_sample_file_n_clicks and not df_g1.empty and not df_g2.empty:
        master = merge_dataframes(df_g1,df_g2)  
        master = protein_create_protein_list(master, species)
        return master
    else:
        return pd.DataFrame()

def create_df_fig(master):
    if not master.empty:
        df_fig = create_protein_df_fig(master)
        return df_fig
    else:
        return pd.DataFrame()

def create_df_info(master, protein_radioitems_value):
    if not protein_radioitems_value:
        protein_radioitems_value == 'area_sum'
    if not master.empty:
        df_protein_info = create_protein_datatable(master, protein_radioitems_value)
        df_protein_info.fillna(0, inplace=True)
        return df_protein_info
    else:
        return pd.DataFrame()


def create_protein_figure_and_table(rows, derived_virtual_selected_rows, search_protein, clickData, protein_radioitems_value, 
checkbox_values, generate_protein_graph_n_clicks, df_fig, df_protein_info, protein_fig):
    if df_fig.empty:
        return {}, start_table_df.to_dict('records'), [{'id': '', 'name': ''}], True, ['Choose protein']
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    highlighted_triv_names = ['Choose protein']
    disabled = True
    if clickData or search_protein or derived_virtual_selected_rows:
        disabled = False
        protein_fig = go.Figure(protein_fig)
        highlighted_triv_names = []
        triv_names_holder = df_fig['trivial_name'].array
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
    if not df_fig.empty:
        if protein_radioitems_value:
            if 'spc_mean' in protein_radioitems_value:
                abundance_metric = 'spc_mean'
                rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','metric_g1': 'SpC G1',  'sd_g1': 'SD 1', 'metric_g2':'SpC G2',  'sd_g2':'SD 2', 'p_val':'p-value'
                }
                
            elif 'area_mean' in protein_radioitems_value:
                abundance_metric = 'area_mean'
                rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','metric_g1': 'intensity G1', 'sd_g1':'SD 1', 'metric_g2':'intensity G2', 'sd_g2': 'SD 2', 'p_val':'p-value'
                }
            elif 'spc_sum' in protein_radioitems_value:
                abundance_metric = 'spc_sum'
                rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2','metric_g1': 'SpC G1',  'sd_g1': 'SD 1', 'metric_g2':'SpC G2',  'sd_g2':'SD 2', 'p_val':'p-value'
                }
            else:
                abundance_metric = 'area_sum'
                rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2', 'metric_g1': 'intensity G1', 'sd_g1':'SD 1', 'metric_g2':'intensity G2',  'sd_g2': 'SD 2', 'p_val':'p-value'
                }
                
        else:
            abundance_metric = 'area_sum'
            rename = {'#peptides_g1': '#peptides G1', '#peptides_g2': '#peptides G2', 'metric_g1': 'intensity G1', 'sd_g1':'SD 1', 'metric_g2':'intensity G2',  'sd_g2': 'SD 2', 'p_val':'p-value'
            }
        
        df_protein_info.rename(columns = rename, inplace=True)
        protein_info_data = df_protein_info.to_dict('rows')
        protein_info_columns=[{"name": str(i), "id": str(i)} for i in df_protein_info.columns]
        if str(changed_id) == 'generate-protein-graph.n_clicks' or str(changed_id) == 'protein-radioitems.value' or str(changed_id) == 'protein-checkbox.value':
            if checkbox_values and 'show-stdev' in checkbox_values :
                protein_fig = create_protein_fig(df_fig, show_stdev = True,  abundance_metric = abundance_metric)
            else:
                protein_fig = create_protein_fig(df_fig, abundance_metric = abundance_metric)
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

def process_peptide_data_for_fig(n_clicks_generate_peptide_fig, peptide_radioitems_value, master, button_label):
    if n_clicks_generate_peptide_fig and not master.empty:
        trivname = button_label.split(' ')[-1]
        peptide_df = master[master['trivname'] == trivname]
        pos_sample, neg_sample, y_label = pre_process_peptide_fig(peptide_df, abundance_metric = peptide_radioitems_value)
        return peptide_df, [pos_sample, neg_sample, trivname, y_label]
    else:
        return pd.DataFrame(), []


def create_peptide_fig_callback(processed_peptide_data, sum_or_mean_radio, rows, derived_virtual_selected_rows):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if changed_id == 'peptide-info-table.derived_virtual_selected_rows' and rows:
        selected_rows_df = pd.DataFrame(rows)
        squares = []
        x0_array, x1_array = list(selected_rows_df.iloc[derived_virtual_selected_rows, 1].values), list(selected_rows_df.iloc[derived_virtual_selected_rows, 2].values)
        for x0, x1 in zip(x0_array, x1_array):
            squares.append((x0, x1))
        if not rows:
            squares = [(0,0)]
    else:
        squares = [(0,0)]
    if processed_peptide_data:
        pos_sample, neg_sample, trivname, y_label = processed_peptide_data
        peptide_fig = create_peptide_fig(pos_sample, neg_sample, trivname, y_label, show_difference='show', show_weight ='show', average=sum_or_mean_radio, square=squares)
        return peptide_fig
    else:
        return {}

def create_peptide_table(peptide_df, peptide_radioitems_value):
    if peptide_radioitems_value == 'area':
        sort = ['intensity G1','intensity G2']
        rename = {'metric_g1':'intensity G1', 'sd_g1':'SD 1','metric_g2':'intensity G2','sd_g2':'SD 2'}
    else:
        sort = ['SpC G1','SpC G2']
        rename = {'metric_g1':'SpC G1', 'sd_g1':'SD 1','metric_g2':'SpC G2','sd_g2':'SD 2'}
    if not peptide_df.empty:
        df_peptide_info = create_peptide_datatable(peptide_df, abundance_metric=peptide_radioitems_value)
        df_peptide_info.fillna(0, inplace=True)
        df_peptide_info.rename(columns=rename, inplace=True)
        df_peptide_info.sort_values(by=sort, ascending=False, inplace=True)
        peptide_table_data = df_peptide_info.to_dict('rows')
        peptide_table_columns=[{"name": str(i), "id": str(i)} for i in df_peptide_info.columns]
        return peptide_table_data, peptide_table_columns
    else: 
        return start_table_df.to_dict('records'), [{'id': '', 'name': ''}]

def create_amino_acid_fig(dropdown_values, radioitem_value, peptide_df, df):
    if dropdown_values and 'complete-proteome' in dropdown_values and not df.empty:
        fig = amino_acid_piecharts(df, accession = '', peptide_or_protein_list = 'protein_list', abundance_metric = radioitem_value)
        return fig
    elif dropdown_values and 'selected-protein' in dropdown_values and not peptide_df.empty:
        accession = peptide_df['Accession'].values[0]
        fig = amino_acid_piecharts(df, accession=accession, peptide_or_protein_list = 'peptide_list', abundance_metric = radioitem_value)
        return fig
    else:
        return  {}


def generate_hover_graphs(hoverData, protein_radioitems_value, master_df):
    if not protein_radioitems_value:
        protein_radioitems_value='area_sum'

    if hoverData:
        accession = hoverData['points'][0]['customdata'][-2]
        fig = all_sample_bar_chart(master_df, accession, metric = protein_radioitems_value)
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

def create_peptide_length_dropdown(length_dropdown_values, amino_acid_radioitems, peptide_df, df):
    if length_dropdown_values and 'complete-proteome'  in length_dropdown_values and not df.empty:
        length_fig = create_length_histogram(df, abundance_metric=amino_acid_radioitems, peptide_or_protein_list='protein_list')
        return length_fig
    elif length_dropdown_values and 'selected-protein' in length_dropdown_values and not peptide_df.empty:
        accession = peptide_df['Accession'].values[0]
        length_fig= create_length_histogram(df, accession=accession, abundance_metric=amino_acid_radioitems, peptide_or_protein_list='peptide_list')
        return length_fig
    else:
        return {}

def create_venn_bar_fig(length_dropdown_values, peptide_df, df):
    if length_dropdown_values and 'complete-proteome'  in length_dropdown_values and not df.empty:
        venn_bar = create_venn_bar(df, accession = '', complete_proteome=True)
        return venn_bar
    elif length_dropdown_values and 'selected-protein' in length_dropdown_values and not peptide_df.empty:
        accession = peptide_df['Accession'].values[0]
        venn_bar= create_venn_bar(df, accession=accession, complete_proteome=False)
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
    State('protein-df-cutoff','data'),
)(generate_hover_graphs)

app.callback(
    Output('aa-fig', 'figure'),
    Input('amino-acid-pie-dropdown', 'value'),
    Input('aa-radioitems','value'),
    State('peptide-df', 'data'),
    State('protein-df-cutoff', 'data'),
)(create_amino_acid_fig)

app.callback(
    Output('peptide-length-fig', 'figure'),
    Input('amino-acid-pie-dropdown', 'value'),
    Input('aa-radioitems','value'),
    State('peptide-df', 'data'),
    State('protein-df-cutoff', 'data'),
)(create_peptide_length_dropdown)

app.callback(
    Output('venn-bar', 'figure'),
    Input('amino-acid-pie-dropdown', 'value'),
    State('peptide-df', 'data'),
    State('protein-df-cutoff', 'data'),
)(create_venn_bar_fig)

app.callback(
    ServersideOutput('peptide-df', 'data'), 
    ServersideOutput('processed-peptide-data', 'data'),
    Input('generate-peptide-fig', 'n_clicks'),
    Input('peptide-radioitems', 'value'),
    State('protein-df-cutoff', 'data'),
    State('generate-peptide-fig', 'children'),
)(process_peptide_data_for_fig)

app.callback(
    Output('peptide-fig', 'figure'),
    Input('processed-peptide-data', 'data'),
    Input('sum-or-mean-radio', 'value'),
    Input('peptide-info-table', "derived_virtual_data"),
    Input('peptide-info-table', 'derived_virtual_selected_rows'),
)(create_peptide_fig_callback)

app.callback(
    Output('peptide-info-table', 'data'),
    Output('peptide-info-table', 'columns'),
    Input('peptide-df', 'data'),
    Input('peptide-radioitems', 'value'),
)(create_peptide_table)

app.callback(
    Output("output-filename-1", "data"),
    Output("sample-collapse-1", "data"),
    ServersideOutput('df_g1-holder', 'data'),
    [Input("upload-data-1", "contents"), Input("upload-data-1", "filename"),
    Input('log-checkbox', 'value')],
)(update_file_list)

app.callback(
    Output("output-filename-2", "data"),
    Output("sample-collapse-2", "data"),
    ServersideOutput('df_g2-holder','data'),
    [Input("upload-data-2", "contents"), Input("upload-data-2", "filename"),
    Input('log-checkbox', 'value')],
)(update_file_list)

app.callback(
    Output('protein-list', 'children'),
    ServersideOutput('protein-df-cutoff', 'data'),
    Input('protein-df', 'data'),
    Input('close-modal-normalization', 'n_clicks'),
    Input('close-modal-cutoff', 'n_clicks'),
    State('cutoff-value-holder', 'children'),
    State('normalization-holder', 'children'),
    State('housekeeping-protein-holder', 'children'),
    )(apply_cutoffs_to_protein_list)

app.callback(
    ServersideOutput('protein-df', 'data'),
    Input("close-modal-file", "n_clicks_timestamp"),
    Input('load-sample-data','n_clicks_timestamp'),
    State('df_g1-holder', 'data'),
    State('df_g2-holder', 'data'),
    State('select-species','value'),
    memoize=True
    )(make_protein_list)


app.callback(
    ServersideOutput('protein-df-fig-holder', 'data'),
    Input('protein-df-cutoff', 'data'),
    memoize=True,
    )(create_df_fig)

app.callback(
    ServersideOutput('protein-datatable-holder','data'),
    Input('protein-df-cutoff', 'data'),
    Input('protein-radioitems','value'),
    memoize=True,
    )(create_df_info)

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
    Input('protein-df-fig-holder', 'data'),
    Input('protein-datatable-holder','data'),
    State('protein-fig', 'figure'),
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
    [Input("open-modal-file", "n_clicks"), Input("close-modal-file", "n_clicks"),
    Input("close-modal-file-2", "n_clicks")],
    [State("modal-file", "is_open")],
)(toggle_modal)


app.callback(
    Output("modal-cutoff", "is_open"),
    [Input("open-modal-cutoff", "n_clicks"), Input("close-modal-cutoff", "n_clicks"),
    Input("close-modal-cutoff-2", "n_clicks")],
    [State("modal-cutoff", "is_open")],
)(toggle_modal)

app.callback(
    Output("modal-normalization", "is_open"),
    [Input("open-modal-normalization", "n_clicks"), Input("close-modal-normalization", "n_clicks"),
    Input("close-modal-normalization-2", "n_clicks")],
    [State("modal-normalization", "is_open")],
)(toggle_modal)


app.callback(
    Output("modal-feedback", "is_open"),
    [Input("open-modal-feedback", "n_clicks"), Input("close-modal-feedback", "n_clicks"),
    Input("close-modal-feedback", "n_clicks")],
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
