import dash_html_components as html
import dash_bootstrap_components as dbc

how_to_use = '''
#### How to use: 


**Step zero:** If you already know what normalizations and cutoffs to apply to your dataset, do so by pressing the **NORMALIZATION** and **CUTOFFS** buttons
 before uploading the files to reduce load time.

**Step one:** Upload your files to the respective groups by pressing the **FILES** button (or load example data by pressing the button below).

**Step two:** Click **GENERATE PROTEIN GRAPH** to inspect the protome. You may click on datapoints in the graph or use the generated table. 
Hovering on a protein will generate a barchart showing the abundance metric of all the samples.

**Step three:** Once a protein is chosen you may generate it's peptide view by clicking the button containing the protein name.

**Step four:** To view the length distribution, amino acid profile, and peptide overlap in your sample; select _Selected protein_ or 
_Complete proteome_ in the dropdown under **General Characteristics** (_selecting complete proteome may result in a long loading time_). 

'''


Data_processing = html.Div([html.H5('Data processing', style={"font-weight": "bold"}),
    html.P('''XXXX requires a minimum of one input file per group, however, statistical calculations require at least 
    three files per group. The input files should either be in a CSV (.csv) or Excel (.xlsx) format. The files are stored 
    in your browser during the session.'''),
    dbc.Card(html.P( '''The input files require one column of: Peptide sequence, UniProt id and intensity measure. 
    { Peptide, Protein, Intensity }'''), color='info', inverse=True),
])

Settings_and_filters = html.Div()

Visualization = html.Div(
)

Interactivity = html.Div()

Cite = html.Div()

Legal = html.Div()

Contact = html.Div()


Documentation = html.Div([
    Data_processing,
    Settings_and_filters,
    Visualization,
    Interactivity,
    Cite,
    Legal,
    Contact,


])