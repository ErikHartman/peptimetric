
import pandas as pd
import plotly.express as px 
import plotly.graph_objects as go

import dash  
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import sys, os

app = dash.Dash(__name__)

# Import files 

# App Layout 

app.layout = html.Div([

    html.H1("Visualizatition tool for proteomics ", style={'text-align': 'center'}), 

    dcc.Textarea(
        id='textarea-example',
        placeholder='Serach protein...', 
        value='Hello World',
        style={'width': '40%'} , 
   ),
    html.Div(id='textarea-example-output', style={'whiteSpace': 'pre-line'})
])

if __name__ == '__main__':
    app.run_server(debug=True)
