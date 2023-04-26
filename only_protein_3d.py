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
from main import adding_cols
pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)

# incorporate data into app
df = pd.read_csv("/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/merged_brca1_sge_ukb_2023_04_21.csv")
exon_list = [(43125364, 43125271), (43124115, 43124017), (43115779, 43115726), (43106533, 43106456),
             (43104956, 43104868), (43104261, 43104122), (43099880, 43099775), (43097289, 43097244),
             (43095922, 43095846), (43094860, 43091435), (43091032, 43090944), (43082575, 43082404),
             (43076614, 43076488), (43074521, 43074331), (43071238, 43070928), (43067695, 43067608),
             (43063951, 43063874), (43063373, 43063333), (43057135, 43057052), (43051117, 43051063),
             (43049194, 43049121), (43047703, 43047643), (43045802, 43044295)]
intron_list = [(43124116, 43125270), (43115780, 43124016), (43106534, 43115725), (43104957, 43106455),
              (43104262, 43104867), (43099881, 43104121), (43097290, 43099774), (43095923, 43097243),
              (43094861, 43095845), (43091033, 43091434), (43082576, 43090943), (43076615, 43082403),
              (43074522, 43076487), (43071239, 43074330), (43067696, 43070927), (43063952, 43067607),
              (43063374, 43063873), (43057136, 43063332), (43051118, 43057051), (43049195, 43051062),
              (43047704, 43049120), (43045803, 43047642)]
df = adding_cols(df, exon_list)

# Build your components------------------------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY],suppress_callback_exceptions=True,
           meta_tags=[{'name': 'viewport',
                            'content': 'width=device-width, initial-scale=1.0'}] )
# 3D parsing & styling
parser = PdbParser('/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/AF-P38398-F1-model_v4.pdb')
# from https://alphafold.ebi.ac.uk/entry/P38398
data = parser.mol3d_data()
styles = create_style_3d(
    df, 'minmax_neg_func_score', data['atoms'], visualization_type='cartoon', color_element='residue_score')
brca1_3D=dashbio.Molecule3dViewer(
                                id='dashbio-default-molecule3d',
                                modelData=data,
                                styles=styles,
                                backgroundOpacity=0,
                                selectionType='residue',
                                backgroundColor="blue",
                                height=900,
                                width=900  # set width to 100%
)

app.layout =brca1_3D


@app.callback(
    Output('default-molecule3d-output', 'children'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    if atom_ids is None or len(atom_ids) == 0:
        phr1='No amino acid has been selected. Click somewhere on the protein \
        structure to select an amino acid.'

        return phr1
    for atm in atom_ids:
        print('Residue name: {}'.format(data['atoms'][atm]['residue_name']))

    aa_name = 'Amino acid: '+data['atoms'][atm]['residue_name']+'\n'
    subset_df=df.loc[df['aa_pos'] == data['atoms'][atm]['residue_index']]
    print(subset_df)
    if len(list(subset_df['var_name'])) < 1:
        return html.Div([html.Br(), html.Div(aa_name),html.Div("No variant correspond to this amino acid"),html.Br()])

    else:
        if len(list(subset_df['var_name'])) == 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variant at this position.' +'\n'
        if len(list(subset_df['var_name'])) > 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variants at this position.'+'\n'
        variant_list=[]
        for i in range(len(list(subset_df['var_name']))):
            variant_list.append(str(list(subset_df['var_name'])[i])+' (SGE function score:'+str(list(subset_df['minmax_neg_func_score'])[i])+')')
        return html.Div([html.Br(), html.Div(aa_name),
            html.Div(variant_nb)]+
            [html.Div(var) for var in variant_list]+
            [html.Br()])


# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8054)
