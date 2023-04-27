# If you prefer to run the code online instead of on your computer click:
# https://github.com/Coding-with-Adam/Dash-by-Plotly#execute-code-in-browser
import numpy as np
from dash import Dash, dcc, Output, Input  # pip install dash
import dash_bootstrap_components as dbc  # pip install dash-bootstrap-components
import plotly.express as px
import pandas as pd  # pip install pandas
from dash import dcc, html, Output, Input
from plotly.subplots import make_subplots

from protein_folding import create_style_3d
import dash_bio as dashbio
from dash import html
from dash_bio.utils import PdbParser
from app import adding_cols
pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)


# add an image

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

image_url = 'https://your-image-url.com/image.jpg'

app.layout = dbc.Container([
    html.Img(src=image_url, style={'height': '300px', 'width': 'auto'})
])