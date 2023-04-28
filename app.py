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

templates = 'plotly_dark'

pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)

from dash.development.base_component import Component, _explicitize_args


# class BooleanSwitch2(Component) copy from https://github.com/plotly/dash-daq/blob/master/dash_daq/BooleanSwitch.py
class BooleanSwitch2(Component):
    """A BooleanSwitch component.
A switch component that toggles
between on and off.
Keyword arguments:
- id (string; optional):
    The ID used to identify this component in Dash callbacks.
- className (string; optional):
    Class to apply to the root component element.
- color (string; optional):
    Color to highlight active switch background.
- disabled (boolean; optional):
    If True, switch cannot be clicked.
- label (dict; optional):
    Description to be displayed alongside the control. To control
    styling, pass an object with label and style properties.
    `label` is a string | dict with keys:
    - label (string; optional)
    - style (dict; optional)
- labelPosition (a value equal to: 'top', 'bottom'; default 'top'):
    Where the component label is positioned.
- on (boolean; default False):
    Whether or not the switch is on.
- persisted_props (list of a value equal to: 'on's; default ['on']):
    Properties whose user interactions will persist after refreshing
    the component or the page. Since only `on` is allowed this prop
    can normally be ignored.
- persistence (boolean | string | number; optional):
    Used to allow user interactions in this component to be persisted
    when the component - or the page - is refreshed. If `persisted` is
    truthy and hasn't changed from its previous value, a `value` that
    the user has changed while using the app will keep that change, as
    long as the new `value` also matches what was given originally.
    Used in conjunction with `persistence_type`.
- persistence_type (a value equal to: 'local', 'session', 'memory'; default 'local'):
    Where persisted user changes will be stored: memory: only kept in
    memory, reset on page refresh. local: window.localStorage, data is
    kept after the browser quit. session: window.sessionStorage, data
    is cleared once the browser quit.
- size (number; optional):
    size of the switch.
- style (dict; optional):
    Style to apply to the root object.
- theme (dict; default light):
    Theme configuration to be set by a ThemeProvider.
- vertical (boolean; default False):
    If True, switch will be vertical instead of horizontal."""

    @_explicitize_args
    def __init__(self, id=Component.UNDEFINED, on=Component.UNDEFINED, color=Component.UNDEFINED,
                 vertical=Component.UNDEFINED, disabled=Component.UNDEFINED, theme=Component.UNDEFINED,
                 label=Component.UNDEFINED, labelPosition=Component.UNDEFINED, className=Component.UNDEFINED,
                 style=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED,
                 persistence_type=Component.UNDEFINED, size=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'className', 'color', 'disabled', 'label', 'labelPosition', 'on', 'persisted_props',
                            'persistence', 'persistence_type', 'size', 'style', 'theme', 'vertical']
        self._type = 'BooleanSwitch'
        self._namespace = 'dash_daq'
        self._valid_wildcard_attributes = []
        self.available_properties = ['id', 'className', 'color', 'disabled', 'label', 'labelPosition', 'on',
                                     'persisted_props', 'persistence', 'persistence_type', 'size', 'style', 'theme',
                                     'vertical']
        self.available_wildcard_properties = []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(BooleanSwitch2, self).__init__(**args)


def adding_cols(df, exons):
    df = df.sort_values("position")
    df = df.rename(columns={"position": "Genomic position", "consequence_ukb": "Consequence", "clinvar": "Clinvar",
                            "func_class_ukb": "UKB", "func_class_sge": "SGE"})
    df = df.replace({"UKB": {"neutral": "Neutral", "lof": "Loss of Function"}})
    df = df.replace({"SGE": {"neutral": "Neutral", "lof": "Loss of Function", "intermediate": "Intermediate"}})
    df = df.replace({"Consequence": {"synonymous": "Synonymous", "missense": "Missense",
                                     "stop_gained": "Stop gained", "splice_region": "Splice region", "intron": "Intron",
                                     "splice_acceptor": "Splice acceptor", "splice_donor": "Splice donor",
                                     "start_lost": "Start lost", "5_utr": "5' UTR"}})
    df["Clinvar"] = df["Clinvar"].replace('absent', 'Absent')
    df['var_index'] = [x for x in range(len(df['Genomic position']))]
    # normalize function score
    df['minmax_neg_func_score'] = [-a for a in df['func_score'].to_list()]
    df['minmax_neg_func_score'] = (df['minmax_neg_func_score'] - df[
        'minmax_neg_func_score'].min()) / (df['minmax_neg_func_score'].max() - df[
        'minmax_neg_func_score'].min())
    df['clinvar_simple'] = "Absent / no clear interpretation"
    df['clinvar_simple'] = df['Clinvar'].map(
        {"Likely pathogenic": "Pathogenic", 'Pathogenic': 'Pathogenic', "Pathogenic/Likely pathogenic": "Pathogenic",
         "Likely benign": "Benign", "Benign": 'Benign'})
    # Get size depend on AC
    df['1/AC'] = [1 / x for x in df['cohort_allele_count'].to_list()]
    df['cadd_score'] = df['cadd_score'] / 50
    # Get cumulative score
    df["cumulative_score"] = (df['1/AC'] + df['minmax_neg_func_score'] + df['cadd_score']) / 3
    print("cumulative=", min(df["cumulative_score"]), max(df["cumulative_score"]))
    # Get Y axis based alt nucleotide var
    df['alt_pos'] = [x for x in df['alt'].map({"A": 1.25, "C": 1, "G": 0.75, "T": 0.5})]
    df['ref_pos'] = [x for x in df['ref'].map({"A": 1.25, "C": 1, "G": 0.75, "T": 0.5})]
    df['var_name'] = "chr" + df['chr'].astype(str) + " " + df['Genomic position'].astype(str) + " " + df['ref'].astype(
        str) + ">" + df['alt'].astype(str)
    df["Exon"] = "Intron"
    for i in range(len(exons)):
        df["Exon"] = np.where((df["Genomic position"] > exons[i][1]) & (df["Genomic position"] < exons[i][0]),
                              "Exon " + str(i + 1), df["Exon"])
    print(df.head())
    print(df.shape)
    return df


font_list = ["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif", "Droid Sans Mono", "Gravitas One",
             "Old Standard TT", "Open Sans", "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
idx_font = 0  # UCL font

# incorporate data into app
df = pd.read_csv(
    "/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/merged_brca1_sge_ukb_2023_04_21.csv")
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
light_gray = 'rgb(140, 130, 121)'  # UCL color
yellow = 'rgb(246, 190, 0)'  # UCL color
yel = "rgb(214,210,196)"  # UCL color
mid_red = 'rgb(147,39,44)'  # UCL color
mid_green = 'rgb(143,153,62)'  # UCL color
mid_purple = 'rgb(80, 7, 120)'  # UCL color
yel_exon = "rgba(214,210,196,0.1)"  # UCL color
dark_gray = 'rgba(41,41,41,1)'  # background

dark_gray_transp = 'rgba(41,41,41,0.85)'
transparent = 'rgba(0,0,0,0) '
exons_color_l1 = ["rgb(86,235,211)", "rgb(106,16,166)", "rgb(97,242,45)", "rgb(194,87,211)", "rgb(56,120,54)",
                  "rgb(147,208,226)", "rgb(51,58,158)", "rgb(189,155,244)", "rgb(33,74,101)", "rgb(55,141,174)",
                  "rgb(191,1,42)", "rgb(239,187,162)", "rgb(120,48,25)", "rgb(239,151,45)"]
exons_color_l2 = ["rgb(33,240,182)", "rgb(127,174,234)", "rgb(179,241,187)", "rgb(239,106,222)", "rgb(127,238,63)",
                  "rgb(250,117,107)", "rgb(65,216,244)", "rgb(189,137,221)", "rgb(203,223,81)", "rgb(144,164,121)",
                  "rgb(246,238,250)", "rgb(200,147,105)", "rgb(254,143,6)", "rgb(243,192,17)"]
# Build your components------------------------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY], suppress_callback_exceptions=True,
           meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}])
# app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP], assets_folder='path/to/assets')

# 3D parsing & styling
parser = PdbParser('/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/AF-P38398-F1-model_v4.pdb')
# from https://alphafold.ebi.ac.uk/entry/P38398
data = parser.mol3d_data()
styles = create_style_3d(
    df, 'minmax_neg_func_score', data['atoms'], visualization_type='cartoon', color_element='residue_score')

brca1_3D = dashbio.Molecule3dViewer(
    id='dashbio-default-molecule3d',
    modelData=data,
    styles=styles,
    backgroundOpacity=0,
    selectionType='residue',
    backgroundColor="black",
    height=500,
    width=1250, zoom=dict(
        factor=1.9, animationDuration=30000, fixedPath=False))
# ------------------------------------------------------------
overview_title = dcc.Markdown(children='', style=dict(font_family=font_list[idx_font], font_color=yel))
overview_display = dcc.RadioItems(options=["Variants aggregated by position", "Variants expanded by nucleotide type"],
                                  value='Variants aggregated by position', labelClassName="custom-text p-3")
items = ['Clinvar', 'Consequence', "SGE", "UKB"]
overview_dropdown = dcc.Dropdown(options=['Clinvar', 'Consequence', "SGE", "UKB"],
                                 value='Consequence', clearable=False, className='my-custom-dropdown')
overview_graph = dcc.Graph(figure={}, config={
    'staticPlot': False,  # True, False
    'scrollZoom': False,  # True, False
    'doubleClick': 'reset',  # 'reset', 'autosize' or 'reset+autosize', False
    'showTips': True,  # True, False
    'displayModeBar': 'hover',  # True, False, 'hover'
    # 'watermark': False,
    'displaylogo': False,
    'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', "zoom2d"]
}, selectedData=None)
color_blind_option = BooleanSwitch2(on=False, size=30,
                                    label=dict(label="color blind friendly", style=dict(font_color=yel)),
                                    color=mid_purple, labelPosition="left")
exon_option = BooleanSwitch2(on=False, size=30, label=dict(label="exon highlighting", style=dict(font_color=yel)),
                             color=yellow, labelPosition="left")
# select Graph
three_d_graph = dcc.Graph(figure={},
                          config={'staticPlot': False, 'scrollZoom': True, 'doubleClick': 'reset', 'showTips': True,
                                  'displayModeBar': False, 'watermark': False})
three_d_title = dcc.Markdown(children='all variant')
#
clinvar_hist_graph_sge = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                      'showTips': True, 'displayModeBar': False, 'watermark': False})
clinvar_hist_graph_ac = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                     'showTips': True, 'displayModeBar': False, 'watermark': False})
clinvar_hist_graph_cadd = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                       'showTips': True, 'displayModeBar': False, 'watermark': False})

mol_viewer_colorbar = dcc.Graph(figure={}, config={'staticPlot': True, 'scrollZoom': False, 'showTips': False,
                                                   'displayModeBar': False, 'watermark': False})

github_link = html.Div([
    html.A(
        id='gh-link',
        children=['View on GitHub'],
        href="http://github.com/",
        style={'color': yel, 'border': yel},
    ),
    html.Img( src='/assets/GitHub-Mark-64px.png',
        style={'height': '50px', 'margin-left': '10px'},
    ),
], style={
    'background': 'black',
    'color': 'white',
    'height': '80px',
    'width': '100%',
    'border': 'solid 2px white',
    'display': 'flex',
    'align-items': 'center',
    'justify-content': 'center',
    'padding': '20px',
})

text_abreviation = dbc.Card(
    [
        dbc.CardImg(src="/assets/ucl-banner-port-stone-rgb-lg.png", top=True),
        dbc.CardBody(
            [
                html.H4("Card title", className="card-title"),
                html.P(
                    "This is the body of the card. You can put any text or HTML content here.",
                    className="card-text",
                ),
                github_link,

                dbc.CardLink("Genome Function Lab", href="https://www.crick.ac.uk/research/labs/greg-findlay/",
                             target="_blank"),
                dbc.CardLink("UKBiobank initiative",
                             href="https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/about-us/", target="_blank"),
            ],className='my-custom-background'
        )
    ],
    style={"height": "500px", 'background-color': 'rgba(41,41,41,0)','border':'solid 2px rgb(214,210,196)', 'border-radius': '20px',
        'overflow': 'hidden'}
)

# Customize your own Layout--------------------------------------------------------------------------------------------------
row_style = {'display': 'flex', 'flex-wrap': 'wrap', 'align-items': 'stretch'}
app.layout = \
    dbc.Container([
        dbc.Row([html.Br()]),

        dbc.Row([
            dbc.Col(html.H1("Variant effect in BRCA1 gene", className='custom-h1'), width={'size': 7, 'offset': 2}, ),
            dbc.Col([
                dbc.Row(color_blind_option, className="my-custom-switch"),
                dbc.Row(exon_option, className="my-custom-switch")], width={'size': 2}, align='right')
        ], justify='between'),
        # row buffer
        dbc.Row([html.Br()]),
        dbc.Row([html.Br()]),

        # Overview graph element ---
        # dbc.Row([
        #     dbc.Col(overview_title, width=12, className='my-custom-title' )
        # ]),
        dbc.Row([dbc.Col(overview_display, width={'size': 6}),
                 dbc.Col([overview_dropdown], width={'size': 2}, align='right', className='my-custom-dropdown'),
                 ], justify='between'),
        dbc.Row([
            dbc.Col(overview_graph, width=12)
        ], justify='around'),

        # row buffer
        dbc.Row([html.Br()]),

        # row 2 ----------------------
        dbc.Row([
            dbc.Col([
                dbc.Card(
                    [
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
                                                            html.H4(className='app-controls-block', children='What is VarEffectViz?'),
                                                            html.P('VarEffectViz is a visualizer that presents variant effect interpretation for the BRCA1 gene using multiple sources of data, including allele count in UKBiobank, the CADD computational score, AlphaFold  protein structure prediction and the results of saturation genome editing experiments. '),
                                                            html.P('By integrating data from these diverse sources, VarEffectViz provides a more comprehensive view of the pathogenicity of each variant. This approach allows researchers, clinicians, and genetics enthusiasts to assess the impact of genetic variation on the function of the BRCA1 gene, which has important implications for cancer risk and prevention.'),
                                                            html.P('The "Glossary" tab explains key concepts related to genetic variants and their effects. This tab is especially useful for non-experts, as it can help them better understand the terminology and concepts used in the visualization board.')
                                                        ])
                                                    ),
                                                    dcc.Tab(
                                                        label='Variant effect',
                                                        value='var-effect',
                                                        children=html.Div(className='control-tab', children=[
                                                            html.H4(className='app-controls-block',
                                                                    children='The variant interpretation challenges'),
                                                            html.P(
                                                                'Every human being has a unique genetic code that is responsible for many of their physical and biological characteristics. The genetic code is made up of DNA, which is organized into distinct units called genes. When changes occur in the DNA sequence, these changes are known as variants. Some variants have no effect on health or function, while others can lead to disease or altered biological processes.'),
                                                            html.P(
                                                                'To better understand the potential effects of a variant, scientists and clinicians use a process called variant interpretation. This involves analyzing the genetic changes to determine whether they are benign or pathogenic. A combination of different approaches is used in this process, including assessing the frequency of the variant in the general population, using computational algorithms to predict the effect of the variant on biological function, and performing experimental studies to validate these predictions.'),
                                                            html.P(
                                                                'However, variant interpretation is not a perfect process, and each approach has its own strengths and limitations. Therefore, to achieve the best possible understanding of genetic variation, researchers and clinicians must combine different approaches to obtain a more comprehensive view of the potential effects of a variant.'),
                                                            html.P(
                                                                'By using multiple sources of data, variant interpretation can provide valuable insights into the role of genetics in health and disease. It can help guide medical diagnosis, treatment, and prevention by identifying genetic changes that may increase the risk of disease or influence response to treatment. Overall, variant interpretation is a crucial tool in the field of genetics and is essential for advancing our understanding of the role of genetics in health and disease.')
                                                        ])
                                                    ),
                                                    dcc.Tab(
                                                        label='Data',
                                                        value='data-resource',
                                                        children=html.Div(className='control-tab', children=[
                                                            html.H4(className='app-controls-block', children='Data sources'),
                                                            html.Li(
                                                                "AlphaFold: A deep learning system developed by Google's DeepMind that predicts the 3D structure of a protein based on its amino acid sequence. The predicted structure can provide valuable information for variant interpretation, as variants can disrupt protein folding and binding, potentially affecting the function of the protein."),
                                                            html.Li(
                                                                "CADD: A tool that predicts the deleteriousness of genetic variants based on their similarity to known pathogenic and benign variants. The score takes into account various genomic annotations, such as conservation, functional genomics, and regulatory information, to provide a single score that reflects the likelihood of a variant being deleterious. The higher the CADD score, the more pathogenic the variant is predicted to be."),
                                                            html.Li(
                                                                "ClinVar: A public database of genetic variants and their clinical significance. ClinVar is considered a highly trusted, \"gold standard\" resource for variant interpretation, as it collects and curates variant data from multiple sources, including clinical laboratories, research studies, and expert panels. However, it's important to note that the majority of genetic variants have not yet been annotated in ClinVar, meaning that they are classified as \"unknown significance\" or \"absent from ClinVar\"."),
                                                            html.Li(
                                                                "Saturation Genome Editing: A lab technique for testing the effects of genetic variants on protein function in cells. It involves systematically introducing mutations into the DNA sequence of a gene using a genome editing tool such as CRISPR-Cas9, and then assessing the resulting changes in the function of the protein that the gene encodes."),
                                                            html.Li(
                                                                "Saturation Genome Editing function score: A score assigned to a genetic variant based on its effect on cell survival rates when the variant is introduced into a cell using the Saturation Genome Editing technique. The score reflects the degree to which the variant disrupts the normal function of the protein encoded by the gene. Saturation Genome Editing function scores are experimental measures of variant pathogenicity."),
                                                            html.Li(
                                                                "UK Biobank: The biggest human sequencing project to date, containing genetic and health-related data from over 500,000 participants in the United Kingdom. UK Biobank is a valuable resource for researchers studying the genetic basis of diseases, as it provides a large sample size and diverse set of genetic and health-related data.")

                                                        ])
                                                    ),

                                                    dcc.Tab(
                                                        label='Glossary',
                                                        value='interpe',
                                                        children=html.Div(className='control-tab', children=[
                                                            html.Tr([html.Td('BRCA1'), html.Td(
                                                                'A gene that encodes a protein involved in DNA repair and maintenance of genomic stability. Mutations in the BRCA1 gene are associated with an increased risk of developing breast and ovarian cancer.')]),
                                                            html.Tr([html.Td('Variant'), html.Td(
                                                                'A genetic variation that occurs in an individual\'s DNA sequence. Variants can be benign or pathogenic and can affect various biological processes in the body.')]),
                                                            html.Tr([html.Td('DNA'), html.Td(
                                                                'The molecule that carries genetic information and is present in almost all living organisms.')]),
                                                            html.Tr([html.Td('Nucleotide'), html.Td(
                                                                'The basic building block of DNA and RNA, consisting of a sugar molecule, a phosphate group, and one of four nitrogenous bases: adenine (A), cytosine (C), guanine (G), or thymine (T) in DNA (or uracil (U) in RNA). The sequence of nucleotides in DNA encodes genetic information, while the sequence of nucleotides in RNA is used to direct protein synthesis.')]),
                                                            html.Tr([html.Td('RNA'), html.Td(
                                                                'A molecule that plays a critical role in the transfer of genetic information from DNA to proteins. RNA is involved in a variety of biological processes, including protein synthesis and regulation of gene expression.')]),
                                                            html.Tr([html.Td('Amino acid'), html.Td(
                                                                'The building blocks of proteins. Amino acids are linked together to form long chains that fold into specific three-dimensional shapes to carry out various biological functions.')]),
                                                            html.Tr([html.Td('Reference genome'), html.Td(
                                                                'A standard DNA sequence used as a reference for comparing genetic variation in different individuals. It provides a basis for identifying genetic differences that may contribute to disease or other traits.')]),
                                                            html.Tr([html.Td('Exon'), html.Td(
                                                                'A coding region of DNA that contains the instructions for making a protein.')]),
                                                            html.Tr([html.Td('Intron'), html.Td(
                                                                'A non-coding section of DNA that separates the coding regions (exons) of a gene. Introns are removed during the process of making RNA from DNA and do not encode proteins, but can play important regulatory roles in gene expression.')]),
                                                            html.Tr([html.Td('Allele count'), html.Td(
                                                                'The number of copies of a particular variant in a given population.')]),
                                                            html.Tr([html.Td('Loss of Function'), html.Td(
                                                                'A genetic variation that results in a partial or complete loss of the normal function of a protein. Such variants are often associated with pathogenicity and can lead to genetic diseases or an increased risk of disease.')]),
                                                            html.Tr([html.Td('5\' UTR'), html.Td(
                                                                'The non-coding region at the beginning of an mRNA molecule, upstream of the start codon, that plays a role in the regulation of gene expression.')]),
                                                            html.Tr([html.Td('Missense'), html.Td(
                                                                'A genetic variation that results in a change in the amino acid sequence of a protein.')]),
                                                            html.Tr([html.Td('Splice acceptor'), html.Td(
                                                                'A region of DNA that signals the end of an exon and is necessary for proper splicing of the pre-mRNA.')]),
                                                            html.Tr([html.Td("Splice donor"),
                                                                                html.Td("A region of DNA that signals the start of an exon and is necessary for proper splicing of the pre-mRNA.")
                                                                            ]),
                                                                            html.Tr([
                                                                                html.Td("Splice region"),
                                                                                html.Td("A region of DNA that is important for the splicing of pre-mRNA.")
                                                                            ]),
                                                                            html.Tr([
                                                                                html.Td("Start lost"),
                                                                                html.Td("A genetic variation that causes the loss of the initiation codon in the coding sequence of a gene.")
                                                                            ]),
                                                                            html.Tr([
                                                                                html.Td("Stop gained"),
                                                                                html.Td("A genetic variation that creates a premature stop codon in the coding sequence of a gene.")
                                                                            ]),
                                                                            html.Tr([
                                                                                html.Td("Synonymous"),
                                                                                html.Td("A genetic variation that does not change the amino acid sequence of a protein.")
                                                        ]),
                                                        ])
                                                    ),
                                                ])
                                            ],
                                            style={'overflow-y': 'auto', 'max-height': '830px'}
                                        )])])],style={'background-color': 'rgba(41,41,41,0)'}, xs=12, sm=12, md=6, lg=3, xl=3),

            dbc.Col([
                dbc.Card([
                    # dbc.CardHeader('Clinvar Histograms'),
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
                    ],className='my-custom-background')
                ])
            ], className='my-custom-background', xs=12, sm=12, md=6, lg=3, xl=3),
            dbc.Col([three_d_graph], xs=12, sm=12, md=12, lg=5, xl=5)
        ], style=row_style, justify='around'),

        # row buffer
        dbc.Row([html.Br()]),

        # row 3 ----------------------
        dbc.Row([
            # 3D protein
            dbc.Col([html.Div([brca1_3D])], xs=12, sm=12, md=6, lg=7, xl=6),
            # color bar scatter plot
            dbc.Col([html.Div([html.Br(), html.Br(), mol_viewer_colorbar, html.Hr(),
                               html.Div(id='default-molecule3d-output', style={'background-color': dark_gray_transp,
                                                                               'position': 'relative',
                                                                               'z-index': '1'
                                                                               })])], className='custom-text_left',
                    xs=12, sm=12, md=3, lg=3, xl=4),
            dbc.Col([text_abreviation], xs=12, sm=4, md=3, lg=2, xl=2),
        ], style=row_style, justify='around'),
    ], fluid=True)


def histogram(x_axis, color_blind):
    if color_blind:
        colors = {'Benign': mid_green, 'Pathogenic': mid_purple}
    else:
        colors = {'Benign': mid_green, 'Pathogenic': mid_red}

    if x_axis == 'minmax_neg_func_score':
        axlis_label = 'SGE fct score'

    elif x_axis == 'cadd_score':
        axlis_label = 'CADD score'
    else:
        axlis_label = x_axis
    df_t = df
    df_t = df_t.rename(columns={'clinvar_simple': "Clinvar high confidence", "minmax_neg_func_score": "SGE fct score",
                                'cadd_score': 'CADD score'})
    fig = px.histogram(df_t, x=axlis_label, color="Clinvar high confidence", color_discrete_map=colors, marginal="rug",
                       hover_data=["var_name", '1/AC', 'cohort_allele_count', 'CADD score', 'SGE fct score'],
                       hover_name="var_name")

    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    fig.update_layout(
        height=260,
        margin=dict(b=0, t=40, l=1, r=2),
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray,
        font_family=font_list[idx_font],
        xaxis=dict(showgrid=False, visible=True, zeroline=True, title=axlis_label),
        yaxis=dict(showgrid=False, visible=True, zeroline=True, title="variant number"),
        font_color=yel)
    if x_axis == '1/AC':
        fig.update_layout(showlegend=True, legend=dict(title='Clinvar:',orientation='h', yanchor='top', y=-0.3, xanchor='left', x=0),)
    else:
        fig.update_layout(showlegend=False)
    fig.update_yaxes(showgrid=False, row=2, col=1)
    fig.update_xaxes(showgrid=False, row=2, col=1)
    fig.update_yaxes(zeroline=False, row=2, col=1)
    fig.update_xaxes(zeroline=False, row=2, col=1)
    if x_axis == "minmax_neg_func_score":
        fig.update_layout(title_font_family=font_list[idx_font],
                          title_text="Pathogenicity classifiers compared",
                          title_y=0.925,
                          title_font_color=yel,
                          title_font_size=18)

    return fig


# Callback allows components to interact--------------------------------------------------------------------------------------------------
@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    # Output(overview_title, 'children'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value'),
    Input(color_blind_option, 'on')
)
def update_overview_graph(column_name, y_axis_nucleotide,
                          color_blind):  # function arguments come from the component property of the Input
    df_temp = df
    c_green_red = [mid_green, '#9bbf85', '#bbd4a6', '#d7dc99', '#fff8c3', '#f8d192', '#f3a66e', '#ec785c', mid_red]
    c_blind_friendly = [mid_green, '#9bbf85', '#bbd4a6', '#e7f7d5', '#fcf2f8', '#f6d3e8', '#d091bb', '#b3589a',
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
    y_axis = [limit[0] + (limit[1] - limit[0]) / 2 for x in df_temp['Genomic position']]
    yaxis_dict = dict(showgrid=False, visible=False)

    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        marker_symb, size_marker, y_axis, height_grph = "square", 8, 'alt_pos', 340
        yaxis_dict = dict(showgrid=False, zeroline= False, title='Nucleotide', tickvals=[0.5, 0.75, 1, 1.25],
                          ticktext=['T', 'G', 'C', 'A'])
    fig = px.scatter(data_frame=df_temp,
                     x=df_temp['Genomic position'],
                     y=y_axis,
                     height=height_grph,
                     color_discrete_map=dict_color_consq,
                     # size="1/AC",
                     color=column_name,
                     custom_data=["SGE",
                                  "UKB", 'Consequence', 'Clinvar',
                                  'cohort_allele_count', 'var_name', 'Exon'],
                     category_orders={'label': list(dict_color_consq.keys())}
                     )

    fig.update_traces(marker=dict(size=size_marker, symbol=marker_symb), selector=dict(mode="markers"),
                      hovertemplate="<br>".join([
                          "<b>%{customdata[5]}</b>",
                          "SGE: %{customdata[0]}",
                          "UKB: %{customdata[1]}",
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
        ref_plot.update_traces(marker=dict(size=size_marker, symbol="square-open", color=yellow),
                               selector=dict(mode="markers"),
                               name='ref', hovertemplate="<br>".join([
                "<b>Reference allele</b>",
                "Position: %{customdata[0]}",
                "Nucleotide: %{customdata[1]}"
            ]))
        fig.add_trace(ref_plot.data[0])

    # add exon
    for i in range(len(exon_list)):
        start, end = exon_list[i][1] - 0.5, exon_list[i][0] + 0.5
        enlarge = 0.4
        texex = "Exon" + str(i + 1)
        pos_xanchor = 'center'
        if i + 1 in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            texex = texex + "   "
        fig.add_shape(type="rect",
                      x0=start, y0=limit[0] - enlarge, x1=end, y1=limit[1] + enlarge,
                      line=dict(color=yel_exon, width=2),
                      fillcolor=yel_exon, layer='below')
        if i + 1 in [23, 10]:
            start = start + 250
        if i + 1 in [17]:
            start = start + 150
        fig.add_annotation(x=start, y=limit[1] + 0.1,
                           text=texex,
                           showarrow=False,
                           font=dict(color=yel, family=font_list[idx_font], size=10), textangle=270,
                           xanchor=pos_xanchor, yanchor='bottom', bgcolor=dark_gray, opacity=1)
    fig.add_shape(type="rect",
                  x0=start, y0=limit[0] - enlarge - 0.2, x1=end, y1=limit[1] + enlarge,
                  line=dict(color=transparent, width=2),
                  fillcolor=transparent, layer='below')
    # add intron
    for i in range(len(intron_list)):
        start, end = intron_list[i][0] - 0.5, intron_list[i][1] + 0.5
        fig.add_shape(type='line', x0=start, x1=end,
                      y0=limit[0] + (limit[1] - limit[0]) / 2, y1=limit[0] + (limit[1] - limit[0]) / 2,
                      line=dict(color=yel_exon, width=2), layer='below')

    fig.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=1, title="Genomic position"),
        yaxis=yaxis_dict,
        font_family=font_list[idx_font],
        legend=dict(orientation='h', yanchor='top', y=1.2, xanchor='left', x=0,
                    title="<b>" + column_name + ' annotation<b>', font_family=font_list[idx_font]),
        font_color=yel)

    return fig.update_layout(
        uirevision=True
    )  # returned objects are assigned to the component property of the Output


@app.callback(
    Output('default-molecule3d-output', 'children'),
    Output(mol_viewer_colorbar, 'figure'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    # Get a color bar
    fig = px.scatter(df, x='rna.score.1', y='rna.score.2', color='minmax_neg_func_score',
                     color_continuous_scale=['#FFFFFF', mid_red],
                     title="Comparison of RNA replicates")

    fig.update_traces(opacity=1)
    fig.update_layout(
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray_transp,
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        yaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        height=300,
        title_font_size=18,
        title_x=0.95,
        font_color=yel,
        font_family=font_list[idx_font],
        xaxis_title="RNA score",
        yaxis_title="RNA score",
        coloraxis_colorbar=dict(
            title="SGE fct score",
            lenmode="pixels", len=200,
            yanchor="top", y=1.5, x=-0.5,
            tickvals=[0, 1],
            tickmode='array',
            ticks="outside",
            ticktext=["Neutral", "LoF"],
            ticklabelposition='outside right'
        ))

    if atom_ids is None or len(atom_ids) == 0:
        phr1 = 'Click somewhere on the protein \
        structure to select an amino acid.'

        return phr1, fig
    for atm in atom_ids:
        print('Residue name: {}'.format(data['atoms'][atm]['residue_name']))

    aa_name = 'Amino acid: ' + data['atoms'][atm]['residue_name'] + '\n'
    subset_df = df.loc[df['aa_pos'] == data['atoms'][atm]['residue_index']]
    print(subset_df)
    if len(list(subset_df['var_name'])) < 1:
        return html.Div(
            [html.Br(), html.Div(aa_name), html.Div("No variant correspond to this amino acid"), html.Br()]), fig

    else:
        if len(list(subset_df['var_name'])) == 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variant at this position.' + '\n'
        if len(list(subset_df['var_name'])) > 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variants at this position.' + '\n'
        variant_list = []
        for i in range(len(list(subset_df['var_name']))):
            variant_list.append(str(list(subset_df['var_name'])[i]) + ' (SGE function score:' + str(
                list(subset_df['minmax_neg_func_score'])[i]) + ')')
        return html.Div([html.Br(), html.Div(aa_name),
                         html.Div(variant_nb)] +
                        [html.Div(var) for var in variant_list] +
                        [html.Br()]), fig


@app.callback(
    Output(component_id=three_d_graph, component_property='figure'),
    Output(component_id=clinvar_hist_graph_sge, component_property='figure'),
    Output(component_id=clinvar_hist_graph_ac, component_property='figure'),
    Output(component_id=clinvar_hist_graph_cadd, component_property='figure'),
    Input(component_id=overview_graph, component_property="selectedData"),
    Input(color_blind_option, 'on'),
    Input(exon_option, 'on')
)
def update_3d_graph(slct_data, color_blind, exon_option):
    fig3 = histogram("minmax_neg_func_score", color_blind)
    fig4 = histogram("cadd_score", color_blind)
    fig5 = histogram("1/AC", color_blind)
    black3dbg = dict(showbackground=True, backgroundcolor=transparent, gridcolor=light_gray, gridwidth=0.5,
                     zeroline=False)

    if exon_option:
        color_exons = 'Exon'
        legend_showing = True
    else:
        color_exons = 'cumulative_score'
        legend_showing = False
    if slct_data is None or slct_data == {'points': []}:
        fig2 = px.scatter_3d(df, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color=color_exons, color_discrete_sequence=exons_color_l1,
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score',
                                          'minmax_neg_func_score'])  # color_continuous_scale=[('rgb(246, 190, 0)'),(),()])
        fig2.update_traces(hovertemplate="<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        if not exon_option:
            fig2.update_traces(marker=dict(size=4, autocolorscale=True, color=df['cumulative_score']))
        else:
            fig2.update_traces(marker=dict(size=4))
        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family=font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_text="Allele frequency, CADD and SGE function score for all variants",
            title_y=0.955,
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region",
                        font=dict(color=yel)),
            showlegend=legend_showing,
            height=813, coloraxis_colorbar=dict(
                title="Cumulative score",
                lenmode="pixels", len=380, thickness=8
            ))

        return fig2, fig3, fig4, fig5

    if slct_data['points'] == []:
        fig2 = px.scatter_3d(title="Please select at least one variant")
        fig2.update_layout(scene=dict(
            bgcolor=dark_gray),
            paper_bgcolor=dark_gray,
            font_family=font_list[idx_font],
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

        dff2 = df[df.var_name.isin(var)]
        fig2 = px.scatter_3d(dff2, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color=color_exons,
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score',
                                          'minmax_neg_func_score'],
                             color_discrete_sequence=exons_color_l1)  # \nin "+str(set(exons)).replace("{", '').replace("'", '').replace("}", ''))

        fig2.update_traces(hovertemplate="<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        if not exon_option:
            fig2.update_traces(marker=dict(size=4, autocolorscale=True, color=df['cumulative_score']))
        else:
            fig2.update_traces(marker=dict(size=4))

        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family=font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_text="Allele frequency, CADD and SGE function score for all variants",
            title_y=0.955,
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region",
                        font=dict(color=yel)),
            showlegend=legend_showing,
            height=785, coloraxis_colorbar=dict(
                title="Cumulative score",
                lenmode="pixels", len=380, thickness=8
            ))

        return fig2.update_layout(uirevision=True), fig3, fig4, fig5


# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8054)
