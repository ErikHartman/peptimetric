import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc

how_to_use = html.Div([
dbc.Row([html.Img(src = './assets/checklist.svg', style={'height':'6%', 'width':'6%'}),
html.H3('How to use', style={'padding-left':5})]),
html.Hr(style={'margin':5}),
dcc.Markdown('''

**0:** If you already know what normalizations and cutoffs to apply to your dataset, do so by pressing the `NORMALIZATION` and `CUTOFFS` buttons
 before uploading the files to reduce load time.

**1:** Upload your files to the respective groups by pressing the `FILES` button (or load example data by pressing the button below).

**2:** Choose a difference metric from the dropdown. Click `GENERATE PROTEIN GRAPH` to inspect the proteome. You may click on datapoints in the graph or use the generated table. 
Hovering on a protein will generate a barchart showing the abundance metric of all the samples.

**3:** Normalize and apply cutoffs to your data by pressing the `NORMALIZATION` and `CUTOFFS` buttons in the navbar. Note that the cutoffs will be applied 
_after_ normalization. 

**4:** Once a protein is chosen you may generate it's peptide view by clicking the button containing the protein name (eg. `HBA_HUMAN`).

**5:** To view the length distribution, amino acid profile, and peptide overlap in your sample; select _Selected protein_ or 
_Complete proteome_ in the dropdown under **General Characteristics** (_selecting complete proteome may result in a long loading time_). 

''')
])

General = dbc.Card(html.P(children=['''
Peptimetric was developed by Erik Hartman and Simon Mahdavi @ Lunds University to help researchers visualize and explore differences in the proteome and peptidome of sample groups.
 The main features of the peptimetric are: normalizing data, applying cutoffs and showcasing the proteome and peptidome of a dataset. There are many interactive elements to allow for easy manipulation and explorations of the users dataset. The webapp
was developed using ''', html.A('Plotly Dash library for Python', href='https://plotly.com/dash/'), ''' and published on the cloud plattoform ''', html.A('Heroku.', href='https://www.heroku.com/')],
 style={'font-weight':'light'}),  color='#F8F8F8', style={'border':0, 'padding-top':10 })

Data_processing = html.Div([dbc.Row([html.Img(src = './assets/computer.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Data processing', style={"font-weight": "bold"}),]),
    html.P(['''Peptimetric requires a minimum of one input file per group, however, statistical calculations require at least 
    three files per group. The input files should be provided in a comma separated values (.csv, recommended) format. Uploading many large
    files may result in a ''', html.A('server timeout error',href ='https://devcenter.heroku.com/articles/error-codes#h12-request-timeout'), '''. If this happens to your data, consider concatenating some of your .csv files for faster uploads. The files are cached and stored 
    in your browser during the session and will be lost after closing or refreshing the tab (''', html.A('Dash data storage).', href='https://dash.plotly.com/dash-core-components/store')]),
    dbc.Card([
        html.P( ['''The input files require four columns describing the ''', html.A('''peptide sequence, precursor protein ID, intensity and spectral count. ''', style={'font-weight':'bold'})
        , '''To accommodate for the use of different dataprocessing softwares we have included the following allowed names (including some capitalization variation) for the various columns:
    '''], style={'padding-top':15}),
    html.P('Peptide sequence: ', style={'font-weight':'bold', 'margin-bottom':0}), html.P('Peptide, Sequence', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Precursor protein ID:', style={'font-weight':'bold','margin-bottom':0}), html.P('Accession, Protein, UniProt id', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Intensity:', style={'font-weight':'bold','margin-bottom':0}), html.P('Area, Intensity', style= {'color':'#c7254e', 'font-family':'monospace'})
    , html.P('Spectral count:', style={'font-weight':'bold','margin-bottom':0}), html.P('Spectral count, SpC, #Feature', style= {'color':'#c7254e', 'font-family':'monospace'})

     ], color='#DFF0D8', style={'border':0, }),
     html.P(['''A local database is used to fetch protein sequences and names. The database consists of a subset of proteomes from ''', html.A('UniProt', href='https://uniprot.org/'),'''. You therefore have to 
     provide a valid UniProt id (accession number) for each peptide in the dataset and select what species to be used in the search when uploading the files. 
     The database was fetched 2021-04 and contains the following species ''', html.I('''Homo sapiens, Sus scrofa (pig), Rattus norvegicus (rat), Cricetulus griseus (hamster),
     Mus musculus (mouse), Danio rerio (zebra fish), Drosophila melanogaster, Caenorhabditis elegans, Candida albicans and Escherichia coli''', className='font-italic'), 
     '''. If you want to analyze data from a species which do not exist in the database, feel free to ''', html.A('contact us ', href='mailto:peptimetric@gmail.com'), ''' and we'll add it as soon as possible!'''], style={'padding-top':15}),
])

Settings = html.Div([dbc.Row([html.Img(src = './assets/settings.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Pre processing data', style={"font-weight": "bold"}),]),
html.H6('Normalization', style={'font-weight':'bold', 'margin-bottom':0, }),
dbc.Row([dbc.Col([html.P('''
Peptimetric accommodates for two ways of normalizing your data: using the global sample intensity, and by using a housekeeping protein. Both methods 
are valid ways of normalizing MS and MSMS data and may reduce the inter-sample biases introduced in sample preparation and loading. 
''', style={}),
html.P('''All the intensity values from the input files are converted to the logarithm (log) scale. If the values present in the input files are 
already in the logarithmic scale the option in the mark the checkbox in the normalizartion popup. The logaritmic scale is generally used when analysing data from MS, as
logaritmic intensities tends to be normally distrubuted.''', style={})], width={'size':7}),
dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/normalization.jpg', style={'height':'100%', 'width':'100%'})), width={'size':4, 'offset':0}),
]),
html.H6('Cutoffs', style={'font-weight':'bold', 'margin-bottom':0, }),
dbc.Row([dbc.Col([html.P('''
MS samples may contain proteins and peptides of very low quality and/or abundancy.
We therefore give you the option to apply cutoffs to your dataset to remove these proteins.
''', style={}),
html.P(' Peptide cutoffs', style={'font-weight':'bold', 'margin-bottom':0, }),
html.P('''Peptide cutoffs allow you to remove peptides with low quality peptides with low intensity and/or low spectral count. We also give you the option
to remove retention time (RT) and/or collission crossection surface (CCS) outliers. Outliers are considered to be situated
three standard deviations from the mean. The peptide cutoffs are applied before the protein cutoffs.
''', style={}),
html.P('Protein cutoffs', style={'font-weight':'bold', 'margin-bottom':0, }),
html.P('''Protein cutoffs allow you to remove proteins with either low total intensity, low total spectral count or a low amount of peptides.
These cuttoffs are applied after the peptide cutoffs are applied. Therefore, if a protein contains many peptides of low quality or abundance
it will be removed from the dataset if the peptide cutoffs are properly applied.
''', style={})], width={'size':7}),
dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/cutoff.jpg', style={'height':'100%', 'width':'100%'})), width={'size':4, 'offset':0}),
]),

])

Visualization = html.Div([
    dbc.Row([html.Img(src = './assets/scatter.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Visualization', style={"font-weight": "bold",})]),
    html.P('''
    Protein View
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P('''
    The protein view gives you an overview of the proteins present in your samples. The size and color indicate the number of peptides found from
    the precursor protein. The axes in the graph correspond to the chosen difference metric for the two groups. You may change the difference metric
    without significant loading time. The standard deviation is off by default and can be toggled. The graph can be manipulated and saved using the modebar. Hovering
    on a datapoint will reveal a tooltip with information about the protein. It will also generate a bar chart showing the chosen difference metric
    for each sample. Selecting a protein highlights it, making it possible to generate a peptide view.
    '''),
    html.P('''The protein table contain information about the number of peptides and chosen difference metric for the 200 proteins with the largets intensity or spectral count. The table is filterable
    and sorted by default on the sum of both groups' difference metric. The table can further be sorted in multiple ways by using the arrows within the column. Selecting a protein highlights the protein in the protein figure. If you select multiple proteins, only the first protein will be selected
    for generating a peptide view. The protein table can be exported as an Excel file (.xlsx) using the "Export" button. 
    '''),
    dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/protein_view1.jpg', style={'height':'100%', 'width':'100%'})), width= {'size':8, 'offset':2}),
    html.P('''
    Peptide View
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),
    html.P('''
    The peptide view showcases the peptides for the selected precursor protein present in the dataset. There are two modes for demonstrating peptide abundancy:
    stacking each sample and viewing them separately or taking the mean of the group. When viewing each sample, the sample may be toggled by clicking 
    on the sample legend. When viewing the mean the standard deviations may be toggled in the same fashion. The figure also contain the weight 
    and difference between the groups which may be toggled on (default) or off.
    ''', style={} ),
    html.P('''
    The peptide table showcase the 200 peptides with the biggest intensity or spectral count difference in the dataset. The table may be sorted and filtered. 
    Selecting a peptide in the table highlights it in the graph. The peptide table can be exported as an Excel file (.xlsx) using the "Export" button. 
    ''', style={} ),

    dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/peptide_view.jpg', style={'height':'100%', 'width':'100%'})), width= {'size':8, 'offset':2}),

    html.P('''
    General Characteristics
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),

    html.P('''
    The general characteristics section contain three figures: length distribution, peptidome overlap and amino acid profile. You may select to view these figures
    for either the complete proteome or the selected protein (visible in the peptide view). Selecting complete proteome may result in long loading times for large datasets.
     As long as a peptide is present in one of the samples in a group they will be included in the visualization of the general characteristics for that group.
    '''),

    dbc.Row([dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/general_info.jpg', style={'height':'100%', 'width':'100%'}))), 
    dbc.Col(dbc.Card(dbc.CardImg(src ='./assets/aa_profile.jpg', style={'height':'100%', 'width':'100%'}))), 
    ])

])



Interactivity = html.Div([
    dbc.Row([html.Img(src = './assets/interact.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Interactivity', style={"font-weight": "bold",})]),
    html.P('''
    Many of the applications and features in peptrimetric are designed to contain some sort of interactivity. Allowing for an dynamic exploration of your dataset and to interact with the graphics
    in an easy manner. 
    ''', style={}), 
    html.P('''
    Search protein
    ''', style={'font-weight':'bold', 'margin-bottom':0, } ),
    html.P(['''
    The search protein function makes it possible to search for all proteins present in the dataset, and those not removed 
    by any cutoffs. A scrollable and clickable list containing the mnemonic protein identification code (trivial name) from ''' , html.A('UniProt', href='https://uniprot.org/'),'''.
     The selected protein will be highlighted in the protein graph. 
    '''], style={} ),
    html.P('''
    Modebars
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),
    html.P(['''
    Peptrimetric uses the built-in ''', html.A('modebars', href='https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar/'),''' from Plotly to allow 
    the user to interact with the graphs. This provides functions such as
    downloading the publication quality plots to .png-format, zooming, autoscale, reseting axes and allows for different hover options.
    '''], style={} ),
    html.P('''
    Peptide view slider 
    ''', style={'font-weight':'bold', 'margin-bottom':0, }),html.P('''
    The protein sequence slider at the bottom of the peptide view allow you to move along and window the protein sequence. 
    ''', style={} ),
])

Cite = html.Div([
    dbc.Row([html.Img(src = './assets/cite.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Cite', style={"font-weight": "bold"}),]),
    dbc.Card(html.P([ '''Please cite the authors: Erik Hartman and Simon Mahdavi, alongside the URL: ''', html.A('peptimetric.herokuapp.com', href='peptimetric.herokuapp.com')], style={'margin-top':'16px'}), color='#F8F8F8', style={'border':0 }),
])

Legal = html.Div([
    dbc.Row([html.Img(src = './assets/legal.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Legal', style={"font-weight": "bold", }),]),
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
    dbc.Row([html.Img(src = './assets/contact.png', style={'height':'4%', 'width':'4%'}),
        html.H5('Contact', style={"font-weight": "bold", }),]),
    contact_text,
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