import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
from dash_html_components.A import A
from dash_html_components.H4 import H4

how_to_use = html.Div([
dbc.Row([html.Img(src = './assets/checklist.jpg'),
html.H3('How to use', style={'padding-left':5})]),
html.Hr(style={'margin':5}),
dcc.Markdown('''

**0:** If you already know what normalizations and cutoffs to apply to your dataset, do so by pressing the `NORMALIZATION` and `CUTOFFS` buttons
 before uploading the files to reduce load time.

**1:** Upload your files .csv to the respective groups by pressing the `FILES` button (or load example data by pressing the button below).

**2:** Choose an abundance metric from the dropdown. Click `GENERATE PROTEIN GRAPH` to inspect the proteome. 
Hovering on a protein will generate a sample graph, displaying the abundance metric of all the samples.

**3:** (If you wish) normalize and apply cutoffs to your data by pressing the `NORMALIZATION` and `CUTOFFS` buttons in the navbar. Note that the cutoffs will be applied 
_after_ normalization. 

**4:** You may use the search bar, table or simply clicking on the dot to select it. 
Once a protein is selected you may generate it's peptide view by clicking the button containing the protein name (eg. `HBA_HUMAN`).

**5:** To view the length distribution, amino acid profile, and peptide overlap in your sample; select _Selected protein_ or 
_Complete proteome_ in the dropdown under **General Characteristics** (_selecting complete proteome may result in a long loading time_). 

''')
])

General = dbc.Card(html.P(children=['''
Peptimetric was developed by Erik Hartman and Simon Mahdavi @ Lunds University to help researchers visualize and explore differences in peptidomic data of sample groups.
 The main features of the peptimetric are: normalizing data, applying cutoffs and interactively showcasing the proteome and peptidome of a dataset. There are many interactive elements to allow for easy manipulation and exploration of the users dataset. The web app
was developed using ''', html.A('Plotly Dash library for Python', href='https://plotly.com/dash/'), ''' and published on the cloud platform ''', html.A('Heroku.', href='https://www.heroku.com/')],
 style={'font-weight':'light'}),  color='#F8F8F8', style={'border':0, 'padding-top':10 })

Data_processing = html.Div([dbc.Row([html.Img(src = './assets/computer.jpg'),
        html.H5('Input files', style={"font-weight": "bold"}),]),
    html.P(['''Peptimetric requires a minimum of one input file per group, however, statistical analyses require at least 
    three files per group. The input files should be provided in a comma separated values (.csv) format. Uploading many large
    files may result in a ''', html.A('server timeout error',href ='https://devcenter.heroku.com/articles/error-codes#h12-request-timeout'), 
    '''. If this happens to your data, consider concatenating some of your .csv files for faster uploads. The files are stored 
    on Heroku's servers during the session and will be lost after closing or refreshing the tab (''', html.A('Dash data storage).', href='https://dash.plotly.com/dash-core-components/store')]),
    dbc.Card([
        html.P( ['''The input files require four columns describing the ''', html.A('''peptide sequence, precursor protein ID, intensity and spectral count. ''', style={'font-weight':'bold'})
        , '''To accommodate for the use of different data processing softwares we have included the following allowed names (including some capitalization variation) for the various columns:
    '''], style={'padding-top':15}),
    html.P('Peptide sequence: ', style={'font-weight':'bold', 'margin-bottom':0}), html.P('Peptide, Sequence', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Precursor protein ID:', style={'font-weight':'bold','margin-bottom':0}), html.P('Accession, Protein, UniProt id', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Intensity:', style={'font-weight':'bold','margin-bottom':0}), html.P('Area, Intensity', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Spectral count:', style={'font-weight':'bold','margin-bottom':0}), html.P('Spectral count, SpC, #Feature', style= {'color':'#c7254e', 'font-family':'monospace'})

     ], color='#DFF0D8', style={'border':0, }),
     html.P(['''A local database is used to fetch FASTA sequences and UniProt mnemonic identifiers (e.g HBA_HUMAN). The database consists of a subset of proteomes from ''', html.A('UniProt', href='https://uniprot.org/'),'''. You therefore have to 
     provide a valid UniProt id (accession number) for each peptide in the dataset and select what species to be used in the search when uploading the files. 
     The database was fetched 2021-04 and contains the following species ''', html.I('''Homo sapiens, Sus scrofa (pig), Rattus norvegicus (rat), Cricetulus griseus (hamster),
     Mus musculus (mouse), Danio rerio (zebra fish), Drosophila melanogaster, Caenorhabditis elegans, Candida albicans and Escherichia coli''', className='font-italic'), 
     '''. If you want to analyze data from a species which do not exist in the database, feel free to ''', html.A('contact us ', href='mailto:peptimetric@gmail.com'), ''' and we'll add it as soon as possible!'''], style={'padding-top':15}),
])

Settings = html.Div([dbc.Row([html.Img(src = './assets/settings.jpg',),
        html.H5('Pre-processing data', style={"font-weight": "bold"}),]),
html.H6('Normalization', style={'font-weight':'bold', 'margin-bottom':0, }),
dbc.Row([dbc.Col([html.P('''
Peptimetric accommodate for two ways of normalizing your data: using the global sample values, and by using a housekeeping protein. Both methods 
are valid ways of normalizing MS and MSMS data and may reduce the inter-sample biases introduced in sample preparation and loading. When normalizing by 
global sample values, each intensity and spectral count value is divided by the sum of all values in the sample and then multiplied by the group mean. When normalizing by housekeeping protein
all values are divided by the sum of the values of the housekeeping protein and then multiplied by the group mean.
''', style={}),
html.P('''All the intensity values from the input files are converted to the 10th logarithm (log10) scale. If the values present in the input files are 
already in the logarithmic scale mark the checkbox in the normalization popup. The logarithmic scale is generally used when analysing intensities from MS, as
logarithmic intensities tend to follow Gaussian (normal) distribution.''', style={})], width={'size':7}),
dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/normalization.jpg', style={'height':'80%', 'width':'80%'})), width={'size':4, 'offset':0}),
]),
html.H6('Cutoffs', style={'font-weight':'bold', 'margin-bottom':0}),
dbc.Row([dbc.Col([html.P('''
MSMS samples may contain proteins and peptides of very low quality and/or abundance.
Peptimetric therefore give you the option to apply cutoffs to your dataset in order to remove these proteins.
''', style={}),
html.P(' Peptide cutoffs', style={'font-weight':'bold', 'margin-bottom':0}),
html.P('''Peptide cutoffs allow you to remove peptides with low intensity and/or low spectral count. We also give you the option
to remove retention time (RT) and/or collision crossection surface (CCS) outliers. Outliers are considered to be situated
three standard deviations from the mean. The peptide cutoffs are applied before the protein cutoffs.
''', style={}),
html.P('Protein cutoffs', style={'font-weight':'bold', 'margin-bottom':0}),
html.P('''Protein cutoffs allow you to remove proteins with either low total intensity, low total spectral count or a low amount of peptides.
These cutoffs are applied after the peptide cutoffs are applied. Therefore, if a protein contains many peptides of low quality or abundance
it will be removed from the dataset (if the peptide cutoffs are properly applied).
''', style={})], width={'size':7}),
dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/cutoff.jpg', style={'height':'80%', 'width':'80%'})), width={'size':4, 'offset':0}),
]),

])

Visualization = html.Div([
    dbc.Row([html.Img(src = './assets/scatter.jpg',),
        html.H5('Visualization', style={"font-weight": "bold"})]),
    html.P('''
    Protein View
    ''', style={'font-weight':'bold', 'margin-bottom':0,'margin-top':10 } ),
    html.P('''
    The protein view gives you an overview of the proteins present in your samples. The size of a dot indicates the number of peptides found from
    the precursor protein. The axes in the graph correspond to the chosen abudance metric for the two respective groups. The color of the dots indicate
    the distance to the diagonal (where the abundance in both groups are equal). You may change the abudance metric
    without significant loading time. The standard deviation is off by default and can be toggled on. 
    '''),
    html.P('''The protein table contains information about the proteins in the samples, including the number of peptides in both groups, as well as the chosen
    abundance metric and standard deviation.
    '''),

    dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/protein_view_fig.png', style={'height':'100%', 'width':'100%'})), width= {'size':8, 'offset':2}),
    html.P('''
    Peptide View
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),

        html.P('''
        The peptide view showcases the peptides for the selected precursor protein present in the dataset. There are two modes for demonstrating peptide abundance:
        stacking each sample and viewing them separately or taking the mean of the group. When viewing each sample, the individual samples may be toggled by clicking 
        on the sample name in the figure legend. When viewing the mean, the standard deviations may be toggled in the same fashion. The figure also contains the weight 
        and difference between the groups which may be toggled on (default) or off.
        ''', style={} ),
        html.P('''
        The peptide table showcases the peptides present in the sample, alongside the start and end position of the sequence in the precursor protein, as well as the sample mean of the
        chosen abundance metric.
        ''', style={} ),
        
    html.P('''
    General Characteristics
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),
    
    html.P('''
    The general characteristics section contains three figures: length distribution, peptidome overlap and amino acid profile. You may select to view these figures
    for either the complete proteome or the selected protein (visible in the peptide view).
    '''),

dbc.Row([
    dbc.Col(dbc.Card([dbc.CardBody(html.H5('Peptide view')), dbc.CardImg(src ='./assets/peptide_view_fig.png', style={'height':'100%', 'width':'100%'})]), width= {'size':4, 'offset':1}),
    dbc.Col(dbc.Card([dbc.CardBody(html.H5('General characteristics')), dbc.CardImg(src ='./assets/general-characteristics.png', style={'height':'100%', 'width':'100%'})]), width= {'size':4, 'offset':2})
])
])



Interactivity = html.Div([
    dbc.Row([html.Img(src = './assets/interact.jpg',),
        html.H5('Interactivity', style={"font-weight": "bold"})]),
    html.P('''
    All features of peptrimetric are designed to be interactive to allow for fast and dynamic exploration of your dataset.
    ''', style={}), 
    html.P('''
    Modebars
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),
    html.P(['''
    Peptrimetric uses the built-in ''', html.A('modebars from Plotly', href='https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar/'),'''  for
    graph interaction. This provides functions such as
    downloading the publication quality plots to .svg-format, zooming, autoscaling and resetting axes.
    '''], style={} ),
    html.P('''
    Selecting a protein
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P(['''
        Proteins may be highlighted in the protein graph by clicking on them, searching for them in the search bar, or selecting them in the protein table. 
        Once a protein is selected, it may be further investigated in the peptide view.
    '''], style={} ),
    html.P('''
    Hovering on a protein
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P(['''
        Hovering on a protein displays a hoverlabel with some information about its metrics. It also generates a bar chart where the abundance metric for each individual sample
        is displayed.
    '''], style={} ),

    html.P('''
    Peptide graph
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P(['''
        All the individual components of the peptide graph may be hidden or shown by clicking on the figure legend. A slider allow you to focus on specific region in the graph.
    '''], style={} ),
    
    html.P('''
    Selecting a peptide
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P(['''
        Peptides may be selected in the peptide table. Once a peptide is selected, the peptide region is highlighted in the peptide graph.
    '''], style={} ),

    
    
])

Cite = html.Div([
    dbc.Row([html.Img(src = './assets/cite.jpg',),
        html.H5('Cite', style={"font-weight": "bold"}),]),
    dbc.Card(html.P([ '''Please cite the creators: Erik Hartman and Simon Mahdavi, alongside the URL: ''', html.A('peptimetric.herokuapp.com', href='peptimetric.herokuapp.com')], style={'margin-top':'16px'}), color='#F8F8F8', style={'border':0 }),
])

Legal = html.Div([
    dbc.Row([html.Img(src = './assets/legal.jpg',),
        html.H5('Legal', style={"font-weight": "bold"}),]),
    html.P(['''
    Peptimetric is free to use for academic purposes and is available at ''', html.A('GitHub', href="https://github.com/ErikHartman/peptimetric"), ''' under an MIT license. 
    No data is saved and no use is monitored using cookies or third party applications. The Dash library is freely available under the MIT license.
    '''])
])


contact_text = html.P(['''
    Have you found any bugs or would you like to file a feature request? ''', html.A('Click here to contact us! ', href='mailto:peptimetric@gmail.com'),
    '''Please describe your request or problem in as great detail as possible. We will get back to you with further questions if necessary.
    '''],)


Contact = html.Div([
    dbc.Row([html.Img(src = './assets/contact.jpg',),
        html.H5('Contact', style={"font-weight": "bold" }),]),
    contact_text,
])

Example_data= html.Div([
    dbc.Row([html.Img(src = './assets/data.jpg',),
        html.H5('Example data', style={"font-weight": "bold" })]),
    html.P(['''By clicking on Load Example Data in the How to Use-window you generate a dataset derived from ''', html.A('Van et al., 2020,', href='https://pubmed.ncbi.nlm.nih.gov/31879271/'), ''' describing peptidomic analysis of urine.
    The study profiles the urinary peptides from youths with type-1-diabetes and the peptidomes of 15 patients with type-1-diabetes (uploaded as group 2) was compared against 15 without (uploaded as group 1). This study is well suited for analysis using Peptimetric as it contains groups of samples and uses a discovery peptidomics approach.''']),
    
    html.P(['''Van JAD, Clotet-Freixas S, Zhou J, Batruch I, Sun C, Glogauer M, Rampoldi L, Elia Y, Mahmud FH, Sochett E, Diamandis EP, Scholey JW, Konvalinka A. Peptidomic Analysis of Urine from Youths with Early Type 1 Diabetes Reveals Novel Bioactivity of Uromodulin Peptides In Vitro. Mol Cell Proteomics. 2020 Mar;19(3):501-517. ''', html.A('doi:10.1074/mcp.RA119.001858', href='https://www.mcponline.org/article/S1535-9476(20)35037-4/fulltext'),'''. Epub 2019 Dec 26. PMID: 31879271; PMCID: PMC7050109.'''])
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
    Example_data,
    html.Hr(),
    Cite,
    html.Hr(),
    Legal,
    html.Hr(),
    Contact,
])