import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc

how_to_use = html.Div([
dbc.Row([html.Img(src = './assets/checklist.png', style={'height':'6%', 'width':'6%'}),
html.H3('How to use', style={'padding-left':5})]),
html.Hr(style={'margin':5}),
dcc.Markdown('''

**0:** If you already know what normalizations and cutoffs to apply to your dataset, do so by pressing the **NORMALIZATION** and **CUTOFFS** buttons
 before uploading the files to reduce load time.

**1:** Upload your files to the respective groups by pressing the **FILES** button (or load example data by pressing the button below).

**2:** Choose a **difference metric** from the dropdown. Click **GENERATE PROTEIN GRAPH** to inspect the proteome. You may click on datapoints in the graph or use the generated table. 
Hovering on a protein will generate a barchart showing the abundance metric of all the samples.

**3:** Once a protein is chosen you may generate it's peptide view by clicking the button containing the protein name.

**4:** To view the length distribution, amino acid profile, and peptide overlap in your sample; select _Selected protein_ or 
_Complete proteome_ in the dropdown under **General Characteristics** (_selecting complete proteome may result in a long loading time_). 

''')
])

General = html.P('''
XXX was developed by Erik Hartman and Simon Mahdavi to help researchers visualize and explore their proteomic and peptidomic data. 
''', style={'padding-left':15, 'padding-right':15, 'font-weight':'light'})

Data_processing = html.Div([dbc.Row([html.Img(src = './assets/computer.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Data processing', style={"font-weight": "bold"}),]),
    html.P('''XXXX requires a minimum of one input file per group, however, statistical calculations require at least 
    three files per group. The input files should either be in a CSV (.csv) or Excel (.xlsx) format. The files are stored 
    in your browser during the session and will be lost after closing or refreshing the tab.''', style={'padding-left':15, 'padding-right':15}),
    dbc.Card([
        html.P( '''The input files require four columns describing the peptide sequence, precursor protein ID, intensity and spectral count. 
        To accommodate for the use of different dataprocessing softwares we have included the following allowed names for the various columns:
    ''', style={'padding-top':15}),
    html.P('Peptide sequence: ', style={'font-weight':'bold', 'margin-bottom':0}), html.P('Peptide, Sequence, sequence, Sequences, sequences, peptide, peptides')
    , html.P('Precursor protein ID:', style={'font-weight':'bold','margin-bottom':0}), html.P('Accession, Protein, protein, accession, uniprot id, UniProt id, Uniprot id')
    , html.P('Intensity:', style={'font-weight':'bold','margin-bottom':0}), html.P('Area, Intensity, area, intensity, intensities')
    , html.P('Spectral count:', style={'font-weight':'bold','margin-bottom':0}), html.P('Spectral count, SPC, SpC, spc, sc, SC, spectral count, #Feature, spectral counts, #Features')

     ], color='#DFF0D8', style={'border':0, 'padding-left':15, 'padding-right':15}),
])

Settings = html.Div([dbc.Row([html.Img(src = './assets/settings.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Settings', style={"font-weight": "bold"}),]),
html.H6('Normalization', style={'font-weight':'bold', 'margin-bottom':0, 'padding-left':15, 'padding-right':15}),
html.P('''
XXX accommodates for two ways of normalizing your data: using the global intensity and by using a housekeeping protein. Both methods 
are valid ways of normalizing MS and MSMS data. 
''', style={'padding-left':15, 'padding-right':15}),
html.H6('Cutoffs', style={'font-weight':'bold', 'margin-bottom':0, 'padding-left':15, 'padding-right':15}),
html.P('''
Often MS samples contain proteins of very low quality and/or abundancy. These oftten contain few peptides, or have very low intensity and spectral
count. We therefore give you the option to apply cutoffs to your dataset. There are two different cutoffs.
''', style={'padding-left':15, 'padding-right':15}), 
html.P(' Peptide cutoffs', style={'font-weight':'bold', 'margin-bottom':0, 'padding-left':15, 'padding-right':15}),
html.P('''Peptide cutoffs allow you to remove peptides with either low intensity or low spectral count. We also give you the option
to remove retention time and/or collission crossection surface outliers. Outliers are considered to be situated
three standard deviations from the mean. The peptide cutoffs are applied before the protein cutoffs.
''', style={'padding-left':15, 'padding-right':15}),
html.P(' Protein cutoffs', style={'font-weight':'bold', 'margin-bottom':0, 'padding-left':15, 'padding-right':15}),
html.P('''Protien cutoffs allow you to remove proteins with either low total intensity, low total spectral count or a low amount of peptides.
These cuttoffs are applied after the peptide cutoffs are applied. Therefore, if a protein contains many peptides of low quality or abundance
it will be removed from the dataset if the peptide cutoffs are properly applied.
''', style={'padding-left':15, 'padding-right':15})

])

Visualization = html.Div([
        html.H5('Visualization', style={"font-weight": "bold",'padding-left':15, 'padding-right':15})],

)

Interactivity = html.Div([
    dbc.Row([html.Img(src = './assets/interact.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Interactivity', style={"font-weight": "bold",'padding-left':15, 'padding-right':15})]),
html.P('''
''')
])

Cite = html.Div([
    dbc.Row([html.Img(src = './assets/cite.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Cite', style={"font-weight": "bold"}),]),
    dbc.Card(html.P( '''Erik Hartman, Simon Mahdavi'''), color='#F8F8F8', style={'border':0, 'padding-left':15, 'padding-right':15}),
html.P('''
''')
])

Legal = html.Div([
    dbc.Row([html.Img(src = './assets/legal.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Legal', style={"font-weight": "bold", 'padding-left':15, 'padding-right':15}),]),
html.P('''
''')
])

Contact = html.Div([
    dbc.Row([html.Img(src = './assets/contact.png', style={'height':'5%', 'width':'5%'}),
        html.H5('Contact', style={"font-weight": "bold", 'padding-left':15, 'padding-right':15}),]),
html.P('''
Have you found any bugs or would you like to file a feature request? Please use the Feedback button, situated in the navbar and
we'll make sure to get back to you.

For other purposes, feel free to contact us directly.
''')
])


Documentation = html.Div([
    General,
    html.Hr(),
    Data_processing,
    html.Hr(),
    Settings,
    html.Hr(),
    Visualization,
    html.Hr(),
    Interactivity,
    html.Hr(),
    Cite,
    html.Hr(),
    Legal,
    html.Hr(),
    Contact,


])