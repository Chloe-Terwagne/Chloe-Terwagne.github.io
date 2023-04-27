# If you prefer to run the code online instead of on your computer click:
# https://github.com/Coding-with-Adam/Dash-by-Plotly#execute-code-in-browser
import numpy as np
from dash import Dash, dcc, Output, Input  # pip install dash
import dash_bootstrap_components as dbc  # pip install dash-bootstrap-components
import plotly.express as px
import pandas as pd  # pip install pandas
from dash import dcc, html, Output, Input
import dash_daq as daq
from protein_folding import create_style_3d
import dash_bio as dashbio
from dash import html
from dash_bio.utils import PdbParser
import base64

import plotly.io as pio
templates='plotly_dark'

pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)


def adding_cols(df, exons):
    df = df.sort_values("position")
    df = df.rename(columns={"position": "Genomic position", "consequence_ukb":"Variant consequence","clinvar": "Clinvar annotation", "func_class_ukb":"UKB function classification", "func_class_sge":"SGE function classification"})
    df = df.replace({"UKB function classification": {"neutral":"Neutral", "lof":"Loss of Function"}})
    df = df.replace({"SGE function classification": {"neutral":"Neutral", "lof":"Loss of Function", "intermediate":"Intermediate"}})
    df = df.replace({"Variant consequence": {"synonymous":"Synonymous", "missense":"Missense",
                                             "stop_gained":"Stop gained", "splice_region":"Splice region", "intron":"Intron", "splice_acceptor":"Splice acceptor", "splice_donor":"Splice donor", "start_lost":"Start lost", "5_utr":"5' UTR"}})
    df["Clinvar annotation"] = df["Clinvar annotation"].replace('absent', 'Absent')
    df['var_index'] = [x for x in range(len(df['Genomic position']))]
    # normalize function score
    df['minmax_neg_func_score'] = [-a for a in df['func_score'].to_list()]
    df['minmax_neg_func_score'] = (df['minmax_neg_func_score'] - df[
        'minmax_neg_func_score'].min()) / (df['minmax_neg_func_score'].max() - df[
        'minmax_neg_func_score'].min())
    df['clinvar_simple'] = "Absent / no clear interpretation"
    df['clinvar_simple']=df['Clinvar annotation'].map({"Likely pathogenic":"Pathogenic", 'Pathogenic':'Pathogenic',"Pathogenic/Likely pathogenic":"Pathogenic", "Likely benign":"Benign","Benign":'Benign'})
    # Get size depend on AC
    df['1/AC'] = [1 / x for x in df['cohort_allele_count'].to_list()]
    # Get Y axis based alt nucleotide var
    df['alt_pos'] = [x for x in df['alt'].map({"A": 1.25, "C": 1, "G": 0.75, "T": 0.5})]
    df['ref_pos'] = [x for x in df['ref'].map({"A": 1.25, "C": 1, "G": 0.75, "T": 0.5})]
    df['var_name'] = "chr" + df['chr'].astype(str) + " " + df['Genomic position'].astype(str) + " " + df['ref'].astype(
        str) + ">" + df['alt'].astype(str)
    df["Exon"]= "Intron"
    for i in range(len(exons)):
        df["Exon"] = np.where((df["Genomic position"]> exons[i][1]) & (df["Genomic position"]<exons[i][0]), "Exon "+str(i+1), df["Exon"])
    print(df.head())
    print(df.shape)
    return df

font_list=["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif", "Droid Sans Mono", "Gravitas One", "Old Standard TT", "Open Sans", "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
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

# Color Model
# color graph overview:  https://www.learnui.design/tools/data-color-picker.html#divergent
light_gray='rgba(205,205,205,1)'
medium_gray='rgba(64,64,64,1)'
dark_gray='rgba(41,41,41,1)' # background
dark_gray_transp = 'rgba(41,41,41,0.85)'
transparent='rgba(0,0,0,0) '

yel="rgb(214,210,196)" # UCL color
mid_red= 'rgb(147,39,44)'# UCL color
mid_green='rgb(143,153,62)'# UCL color
mid_purple = 'rgb(80, 7, 120)'# UCL color
yel_exon="rgba(214,210,196,0.1)" # UCL color
exons_color_l1= ["rgb(86,235,211)", "rgb(106,16,166)", "rgb(97,242,45)", "rgb(194,87,211)", "rgb(56,120,54)", "rgb(147,208,226)", "rgb(51,58,158)", "rgb(189,155,244)", "rgb(33,74,101)", "rgb(55,141,174)", "rgb(191,1,42)", "rgb(239,187,162)", "rgb(120,48,25)", "rgb(239,151,45)"]
exons_color_l2=["rgb(33,240,182)", "rgb(127,174,234)", "rgb(179,241,187)", "rgb(239,106,222)", "rgb(127,238,63)", "rgb(250,117,107)", "rgb(65,216,244)", "rgb(189,137,221)", "rgb(203,223,81)", "rgb(144,164,121)", "rgb(246,238,250)", "rgb(200,147,105)", "rgb(254,143,6)", "rgb(243,192,17)"]
#exons_color_l1=["rgb(65,187,197)", "rgb(95,112,204)", "rgb(253,146,250)", "rgb(198,62,171)", "rgb(195,177,221)", "rgb(81,131,66)", "rgb(30,212,107)", "rgb(253,44,59)", "rgb(157,199,33)", "rgb(185,92,61)", "rgb(252,162,131)", "rgb(134,116,102)", "rgb(242,172,24)", "rgb(149,114,6)"]
# Build your components------------------------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY],suppress_callback_exceptions=True,
           meta_tags=[{'name': 'viewport','content': 'width=device-width, initial-scale=1.0'}] )
# app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP], assets_folder='path/to/assets')

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
                                backgroundColor="black",
                                height=500,
                                width=1250, zoom=dict(
        factor= 1.9, animationDuration= 30000, fixedPath=False))
idx_font=6
#------------------------------------------------------------
overview_title = dcc.Markdown(children='', style=dict(font_family= font_list[idx_font],font_color=yel))
overview_display = dcc.RadioItems(options=["Variants aggregated by position", "Variants expanded by nucleotide type"],  value='Variants aggregated by position', labelClassName="custom-text p-3" )
items = ['Clinvar annotation', 'Variant consequence', "SGE function classification", "UKB function classification"]
overview_dropdown = dcc.Dropdown(options=['Clinvar annotation', 'Variant consequence', "SGE function classification", "UKB function classification"],
                                value='Variant consequence', clearable=False, className='my-custom-dropdown')
overview_graph = dcc.Graph(figure={}, config={
                      'staticPlot': False,     # True, False
                      'scrollZoom': False,      # True, False
                      'doubleClick': 'reset',  # 'reset', 'autosize' or 'reset+autosize', False
                      'showTips': True,       # True, False
                      'displayModeBar': 'hover',  # True, False, 'hover'
                      #'watermark': False,
                      'displaylogo': False,
                      'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', "zoom2d"]
                        }, selectedData=None)
color_blind_option = daq.BooleanSwitch(on=False, label=dict(label="Color blind friendly", style=dict(font_color=yel)),
                                       color=mid_purple, labelPosition="right")

# select Graph
three_d_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': True,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
three_d_title = dcc.Markdown(children='all variant')
#
clinvar_hist_graph_sge = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
clinvar_hist_graph_ac = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
clinvar_hist_graph_cadd = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})

mol_viewer_colorbar = dcc.Graph(figure={}, config={'staticPlot': True, 'scrollZoom': False,'showTips': False,'displayModeBar': False,'watermark': False})

github_link=html.Div([html.A(id='gh-link',
                        children=[
                            'View on GitHub'
                        ],
                        href="http://github.com/",
                        style={'color': yel, 'border': yel}),
                  html.Img(
                      src='data:image/png;base64,{}'.format(
                          base64.b64encode(
                              open(
                                  '/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/assets/GitHub-Mark-64px.png',
                                  'rb'
                              ).read()
                          ).decode()
                      ))],
                 style={
                     'background': dark_gray,
                     'color': yel,
                 })

text_abreviation = dbc.Card(
    [
        dbc.CardBody(
            [
                html.H4("Card title", className="card-title"),
                html.P(
                    "This is the body of the card. You can put any text or HTML content here.",
                    className="card-text",
                ),
                github_link,
                dbc.CardLink("Genome Function Lab", href="https://www.crick.ac.uk/research/labs/greg-findlay/", target="_blank"),
                dbc.CardLink("UKBiobank initiative", href="https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/about-us/", target="_blank"),
            ]
        )
    ],
    style={"height": "500px"}
)



# Customize your own Layout--------------------------------------------------------------------------------------------------
row_style = {'display': 'flex', 'flex-wrap': 'wrap', 'align-items': 'stretch'}
app.layout = \
    dbc.Container([
        dbc.Row([html.Br()]),

        dbc.Row(
            dbc.Col(html.H1("Variant effect interpretation in BRCA1 gene", className='custom-h1'),
                    width=12)
        ),
        # row buffer
        dbc.Row([html.Br()]),
        dbc.Row([html.Br()]),

        # Overview graph element ---
        dbc.Row([
            dbc.Col(overview_title, width=12, className='my-custom-title' )
        ]),
        dbc.Row([dbc.Col(overview_display, width=7),
                    dbc.Col(overview_dropdown, width=2, className='my-custom-dropdown'),
                    dbc.Col(color_blind_option, width=2, className="my-custom-switch"),
                    ], justify='between'),
        dbc.Row([
                 dbc.Col(overview_graph, width=12)
                 ], justify='around'),

        # row buffer
        dbc.Row([html.Br()]),

        # row 2 ----------------------
        dbc.Row([
            dbc.Col([
                html.Div(
                    id='body',
                    className='app-body',
                    children=[
                        html.Div(
                            id='desc-control-tabs',
                            className='control-tabs',
                            children=[
                                dcc.Tabs(id='about-tabs', value='what-is', children=[
                                    dcc.Tab(
                                        label='About',
                                        value='what-is',
                                        children=html.Div(className='control-tab', children=[
                                            html.H4(className='what-is', children='What is Molecule3D?'),
                                            html.P('Molecule3D is a visualizer that allows you '
                                                   'to view biomolecules in multiple representations: '
                                                   'sticks, spheres, and cartoons.'),
                                            html.P('You can select a preloaded structure, or upload your own, '
                                                   'in the "Data" tab. A sample structure is also '
                                                   'available to download.'),
                                            html.P('In the "View" tab, you can change the style and '
                                                   'coloring of the various components of your molecule.')
                                        ])
                                    ),
                                        dcc.Tab(
                                            label='Data',
                                            value='upload-download',
                                            children=html.Div(className='control-tab', children=[
                                                html.H4(className='app-controls-block', children='How the Data come from'),
                                                html.P('Molecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '                                               
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '                                             
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '                                             
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in muMolecule3D is a visualizer that allows you '
                                                       'to view biomolecules in multiple representations: '
                                                       'sticks, spheres, and cartoons.'),
                                                html.P('You can select a preloaded structure, or upload your own, '
                                                       'in the "Data" tab. A sample structure is also '
                                                       'available to download.'),
                                                html.P('In the "View" tab, you can change the style and '
                                                       'coloring of the various components of your molecule.')
                                            ])
                                    ),
                                    dcc.Tab(
                                        label='inter',
                                        value='interpe',
                                        children=html.Div(className='control-tab', children=[
                                            html.H4(className='what-is', children='What is Molecule3D?'),
                                            html.P('Molecule3D is a visualizer that allows you '
                                                   'to view biomolecules in multiple representations: '
                                                   'sticks, spheres, and cartoons.'),
                                            html.P('You can select a preloaded structure, or upload your own, '
                                                   'in the "Data" tab. A sample structure is also '
                                                   'available to download.'),
                                            html.P('In the "View" tab, you can change the style and '
                                                   'coloring of the various components of your molecule.')
                                        ])
                                    ),
                                ])
                            ],
                            style={'overflow-y': 'auto', 'max-height': '750px'}
                        )])], xs=12, sm=12, md=6, lg=3, xl=3),
            dbc.Col([
                dbc.Card([
                    #dbc.CardHeader('Clinvar Histograms'),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_sge),
                        ]),
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_ac),
                        ]),
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_cadd),
                        ])
                    ])
                ])
            ], xs=12, sm=12, md=6, lg=3, xl=3),
            dbc.Col([three_d_graph], xs=12, sm=12, md=12, lg=5, xl=5)
        ],style=row_style, justify='around'),

        # row buffer
        dbc.Row([html.Br()]),

        # row 3 ----------------------
        dbc.Row([
            #3D protein
            dbc.Col([html.Div([brca1_3D])],  xs=12, sm=12, md=6, lg=7, xl=6),
            #color bar scatter plot
            dbc.Col([html.Div([html.Br(),html.Br(),mol_viewer_colorbar,html.Hr(),html.Div(id='default-molecule3d-output',  style={'background-color': dark_gray_transp,
                                                                                                                        'position': 'relative',
                                                                                                                        'z-index': '1'
                                                                                                                        })])], className='custom-text_left', xs=12, sm=12, md=3, lg=3, xl=4),
            dbc.Col([text_abreviation],  xs=12, sm=4, md=3, lg=2, xl=2),
        ],style=row_style, justify='around'),
    ], fluid=True)


def histogram(x_axis):
    if x_axis=='minmax_neg_func_score':
        axlis_label = 'SGE fct score'
    elif x_axis=='cadd_score':
        axlis_label = 'CADD score'
    else:
        axlis_label=x_axis
    df_t=df
    df_t = df_t.rename(columns={'clinvar_simple':"Clinvar high confidence", "minmax_neg_func_score":"SGE fct score", 'cadd_score':'CADD score' })
    fig = px.histogram(df_t, x=axlis_label, color="Clinvar high confidence", color_discrete_map={'Benign':mid_green, 'Pathogenic':mid_red}, marginal="rug", hover_data = ["var_name", '1/AC', 'cohort_allele_count', 'CADD score', 'SGE fct score'], hover_name="var_name")
    # fig.update_traces(hovertemplate="<br>".join([
    #     "<b>%{hover_data[0]}</b>",
    #     "Exon: %{hover_data[1]}",
    #     "SGE function score: %{hover_data[4]}",
    #     "CADD score: %{hover_data[3]}",
    #     "Number of allele count in UKB: %{hover_data[2]}",
    # ]))

    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    fig.update_layout(
        height=260,
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray,
        font_family=font_list[idx_font],
        xaxis=dict(showgrid=False, visible=True,  zeroline=True, title=axlis_label),
        yaxis=dict(showgrid=False, visible=True, zeroline=True, title="variant number"),
        font_color=yel)
    if x_axis=='1/AC':
        fig.update_layout(showlegend=True, legend=dict(title='Clinvar'))
    else:
        fig.update_layout(showlegend=False)
    fig.update_yaxes(showgrid=False, row=2, col=1)
    fig.update_xaxes(showgrid=False, row=2, col=1)
    fig.update_yaxes(zeroline=False, row=2, col=1)
    fig.update_xaxes(zeroline=False, row=2, col=1)

    return fig

# Callback allows components to interact--------------------------------------------------------------------------------------------------
@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    Output(overview_title, 'children'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value'),
    Input(color_blind_option, 'on')
)
def update_overview_graph(column_name, y_axis_nucleotide, color_blind):  # function arguments come from the component property of the Input
    df_temp = df
    c_green_red = ['#488f31', '#7da84f', '#acc272', '#d7dc99', '#fff8c3', '#f8d192', '#f3a66e', '#ec785c', '#de425b']
    c_blind_friendly = ['#71915e', '#9bbf85', '#bbd4a6', '#e7f7d5', '#fcf2f8', '#f6d3e8', '#d091bb', '#b3589a',
                        mid_purple]
    if color_blind:
        colors = c_blind_friendly
    else:
        colors = c_green_red
    dict_color_consq = {"Neutral": colors[0], "Intermediate": colors[4], "Loss of Function": colors[-1],
                        "Synonymous": colors[0], "Intron": colors[1], "5' UTR": colors[2],
                        "Splice region": colors[3], "Missense": colors[4], "Splice acceptor": colors[5],
                        "Splice donor": colors[6], "Start lost": colors[7], "Stop gained": colors[8],
                        "Benign": colors[0], "Likely benign": colors[1],
                        "Uncertain significance": colors[3], "Absent": colors[4],
                        "Conflicting interpretations of pathogenicity": colors[5],
                        "Likely pathogenic": colors[6], "Pathogenic/Likely pathogenic": colors[7],
                        "Pathogenic": colors[8]}

    list_temp = [x for x in list(dict_color_consq.keys()) if x in list(df[column_name])]

    # Create a categorical variable with the desired order
    df_temp[column_name] = pd.Categorical(df[column_name], categories=list_temp, ordered=True)
    # Sort the dataframe based on the categorical variable
    df_temp = df_temp.sort_values(column_name)

    size_marker, height_grph, marker_symb = 8, 340, "square"
    limit = (0.5, 1.25)
    y_axis = [limit[0] + (limit[1]- limit[0])/2 for x in df_temp['Genomic position']]
    yaxis_dict = dict(showgrid=False, visible=False)

    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        marker_symb, size_marker, y_axis, height_grph = "square", 8, 'alt_pos', 340
        yaxis_dict = dict(showgrid=False, visible=True, title='Nucleotide', tickvals=[0.5, 0.75, 1, 1.25],
                          ticktext=['T', 'G', 'C', 'A'])
    fig = px.scatter(data_frame=df_temp,
                     x=df_temp['Genomic position'],
                     y=y_axis,
                     height=height_grph,
                     color_discrete_map=dict_color_consq,
                     # size="1/AC",
                     color=column_name,
                     custom_data=["SGE function classification",
                     "UKB function classification", 'Variant consequence', 'Clinvar annotation',
                                  'cohort_allele_count', 'var_name', 'Exon'],
                     category_orders={'label': list(dict_color_consq.keys())}
                     )

    fig.update_traces(marker=dict(size=size_marker, symbol=marker_symb), selector=dict(mode="markers"),
                      hovertemplate="<br>".join([
                          "<b>%{customdata[5]}</b>",
                          "SGE classification: %{customdata[0]}",
                          "UKB classification: %{customdata[1]}",
                          "Clinvar classication: %{customdata[3]}",
                          "Consequence: %{customdata[2]}",
                          "Number of allele count in UKB: %{customdata[4]}",
                      ]))

    # add reference variant
    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        ref_plot = px.scatter(data_frame=df_temp,
                              x=df_temp['Genomic position'],
                              y=df_temp['ref_pos'],
                              height=height_grph,
                              custom_data=["Genomic position", "ref", "Exon"],
                              category_orders={'label': list(dict_color_consq.keys())}
                              )
        ref_plot.update_traces(marker=dict(size=size_marker, symbol="square-open"), selector=dict(mode="markers"),
                               name='ref', hovertemplate="<br>".join([
                          "<b>Reference allele</b>",
                          "Position: %{customdata[0]}",
                            "Nucleotide: %{customdata[1]}"
            ]))
        fig.add_trace(ref_plot.data[0])

    # add exon
    for i in range(len(exon_list)):
        start, end = exon_list[i][1] - 0.5, exon_list[i][0] + 0.5
        enlarge=0.4
        texex="Exon" + str(i + 1)
        pos_xanchor = 'center'
        if i + 1 in [1,2,3,4,5,6,7,8,9]:
            texex = texex+"   "
        fig.add_shape(type="rect",
                      x0=start, y0=limit[0]-enlarge, x1=end, y1= limit[1]+enlarge,
                      line=dict(color=yel_exon, width=2),
                      fillcolor=yel_exon, layer='below')
        if i + 1 in [23,10]:
            start = start +250
        if i + 1 in [ 17]:
            start = start +150
        fig.add_annotation(x=start, y=limit[1]+0.1,
                           text=texex,
                           showarrow=False,
                           font=dict(color=yel,family=font_list[idx_font], size=8), textangle=270, xanchor=pos_xanchor, yanchor='bottom', bgcolor=dark_gray, opacity=1)

    # add intron
    for i in range(len(intron_list)):
        start, end = intron_list[i][0] - 0.5, intron_list[i][1] + 0.5
        fig.add_shape(type='line', x0=start, x1=end,
                      y0=limit[0] + (limit[1] - limit[0]) / 2, y1=limit[0] + (limit[1] - limit[0]) / 2,
                      line=dict(color=yel_exon, width=2), layer='below')

    fig.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        xaxis=dict(showgrid=False, visible=False, linecolor=None, linewidth=1, title="Genomic position"),
        yaxis=yaxis_dict,
        font_family=font_list[idx_font],
        legend=dict(orientation='h', yanchor='top', y=1.5, xanchor='left', x=0, title="Annotation", font_family=font_list[idx_font]),
        font_color = yel)


    return fig.update_layout(
    uirevision=True
), '#### BRCA1 gene annotated with ' + column_name.lower()  # returned objects are assigned to the component property of the Output


# define a white to red gradient
#colors = ['#FFFFFF', '#FF0000']

@app.callback(
    Output('default-molecule3d-output', 'children'),
    Output(mol_viewer_colorbar, 'figure'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    # Get a color bar
    fig=px.scatter(df, x='rna.score.1', y='rna.score.2', color='minmax_neg_func_score', color_continuous_scale=['#FFFFFF', mid_red],
                   title="Comparison of RNA replicates")

    fig.update_traces(opacity=1)
    fig.update_layout(
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray_transp,
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        yaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        height=300,
        title_font_size= 18,
        title_x=0.95,
        font_color=yel,
        font_family=font_list[idx_font],
        xaxis_title="RNA score",
        yaxis_title="RNA score",
        coloraxis_colorbar=dict(
            title="SGE fct score",
            # thicknessmode="pixels", thickness=50,
            lenmode="pixels", len=200,
            yanchor="top", y=1.5,x=-0.5,
            tickvals=[0,1],
            tickmode='array',
            ticks="outside",
            ticktext=["Neutral", "LoF"],
            #ticks="outside",
            ticklabelposition='outside right'

            # dtick=5
        ))

    if atom_ids is None or len(atom_ids) == 0:
        phr1='Click somewhere on the protein \
        structure to select an amino acid.'

        return phr1, fig
    for atm in atom_ids:
        print('Residue name: {}'.format(data['atoms'][atm]['residue_name']))

    aa_name = 'Amino acid: '+data['atoms'][atm]['residue_name']+'\n'
    subset_df=df.loc[df['aa_pos'] == data['atoms'][atm]['residue_index']]
    print(subset_df)
    if len(list(subset_df['var_name'])) < 1:
        return html.Div([html.Br(), html.Div(aa_name),html.Div("No variant correspond to this amino acid"),html.Br()]), fig

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
            [html.Br()]), fig

@app.callback(
    Output(component_id=three_d_graph, component_property='figure'),
    Output(component_id=clinvar_hist_graph_sge, component_property='figure'),
    Output(component_id=clinvar_hist_graph_ac, component_property='figure'),
    Output(component_id=clinvar_hist_graph_cadd, component_property='figure'),
    Input(component_id=overview_graph, component_property="selectedData")
)



def update_3d_graph(slct_data):
    fig3 = histogram("minmax_neg_func_score")
    fig4 = histogram("cadd_score")
    fig5 = histogram("1/AC")
    black3dbg = dict(
        showbackground=True,
        backgroundcolor=transparent,
        gridcolor=light_gray,
        gridwidth=0.5,
        zeroline=False)

    if slct_data is None or slct_data == {'points': []}:
        fig2 = px.scatter_3d(df, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color='Exon',
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score', 'minmax_neg_func_score'],
                             color_discrete_sequence=exons_color_l1)
        fig2.update_traces(hovertemplate = "<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family= font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_text= "Allele frequency, CADD and SGE function score for all variants",
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region", font=dict(color=yel )),
            height=785)

        return fig2, fig3, fig4, fig5

    if slct_data['points'] == []:
        fig2 = px.scatter_3d(title= "Please select at least one variant")
        fig2.update_layout(scene=dict(
            bgcolor=dark_gray),
            paper_bgcolor=dark_gray,
            font_family= font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_font_color=yel,
            title_font_size=18,
            height=785)
        return fig2, fig3, fig4, fig5
    else:
        print(f'selected data: {slct_data}')
        exons = [slct_data['points'][i]['customdata'][-1] for i in range(len(slct_data['points']))]
        var = [slct_data['points'][i]['customdata'][-2] for i in range(len(slct_data['points']))]

        print(set(exons))
        dff2 = df[df.var_name.isin(var)]
        fig2 = px.scatter_3d(dff2, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color='Exon',
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score', 'minmax_neg_func_score'],
                             color_discrete_sequence=exons_color_l1, title="Allele frequency, CADD and SGE function score for a subset of variants")# \nin "+str(set(exons)).replace("{", '').replace("'", '').replace("}", ''))
        fig2.update_traces(hovertemplate = "<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family= font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region", font=dict(color=yel )),
            height=785)


        return fig2, fig3, fig4, fig5





# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8054)
