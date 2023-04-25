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
import matplotlib.colors as mcolors

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
app = Dash(__name__, external_stylesheets=[dbc.themes.LUX],suppress_callback_exceptions=True)

# 3D parsing & styling
parser = PdbParser('/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/AF-P38398-F1-model_v4.pdb')
# from https://alphafold.ebi.ac.uk/entry/P38398
data = parser.mol3d_data()
styles = create_style_3d(
    df, 'minmax_neg_func_score', data['atoms'], visualization_type='cartoon', color_element='residue_score')
#------------------------------------------------------------
overview_title = dcc.Markdown(children='')
overview_display = dcc.RadioItems(options=["Variants aggregated by position", "Variants expanded by nucleotide type"], value='Variants aggregated by position')
overview_dropdown = dcc.Dropdown(options=['Clinvar annotation', 'Variant consequence', "SGE function classification", "UKB function classification"],
                                 value='Variant consequence',  # initial value displayed when page first loads
                                 clearable=False)
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

# select Graph
three_d_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': True,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
three_d_title = dcc.Markdown(children='all variant')
#
clinvar_hist_graph_sge = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
clinvar_hist_graph_ac = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})
clinvar_hist_graph_cadd = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': False,'watermark': False})

mol_viewer_colorbar = dcc.Graph(figure={}, config={'staticPlot': True, 'scrollZoom': False,'showTips': False,'displayModeBar': False,'watermark': False})
text_area =dcc.Tabs(id='my-tabs', value='tab-1', children=[
        dcc.Tab(label='About', value='tab-1', children=[
            dcc.Textarea(
                value='This is the About tab.',
                contentEditable = False,
                style={'width': '100%', 'height': '100%', 'resize': 'none'}
            )
        ]),
        dcc.Tab(label='Glossary', value='tab-2', children=[
            dcc.Textarea(
                value='This is the Glossary tab.',
                style={'width': '100%', 'height': '300px'},
                contentEditable=False
            )
        ]),
    dcc.Tabs(id='mol3d-tabs', value='what-is', children=[
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

    ]),
        dcc.Tab(label='Other', value='tab-4', children=[
            dcc.Textarea(
                value='This is the Other tab.',
                style={'width': '100%', 'height': '300px'},
                contentEditable=False
            )
        ])
    ])
text_abreviation= dcc.Tab(
            label='Legen',
            value='what-is',
            children=html.Div(className='control-tab', children=[
                html.H4(className='what-is', children='What is Molecule3D?'),
                html.P('Molecule3D is a visualizer that allows you '
                       'to view biomolecules in multiple representations: '
                       'sticks, spheres, and cartoons.')
            ]))
# Customize your own Layout--------------------------------------------------------------------------------------------------
app.layout = html.Div([
    html.H1("Exploration of BRCA1", style={'text-align': 'center'}),
    html.Label('Grid Layout'),
        dbc.Row([dbc.Col([
                            dbc.Card([
                                dbc.CardHeader(overview_title),
                                dbc.CardBody([
                                    html.Div([
                                        dbc.Row([
                                            dbc.Col(overview_display, width=6),
                                            dbc.Col(overview_dropdown, width=6)
                                        ]),
                                        dbc.Row([
                                            dbc.Col(overview_graph, width=12)
                                        ])
                                    ])
                                ])
                            ])
                        ], width=12)
                    ]),
        dbc.Row([
            dbc.Col(text_area, width=3),
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader('Clinvar Histograms'),
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
            ], width=3),

            dbc.Col(three_d_graph, width=6)
        ]),
        dbc.Row([
        #dbc.Col(dbc.Alert('Another component', color='info'), width=3),
        dbc.Col(html.Div([
                        dashbio.Molecule3dViewer(
                            id='dashbio-default-molecule3d',
                            modelData=data,
                            styles=styles,
                            backgroundOpacity=0,
                            selectionType='residue',
                            backgroundColor="blue",
                         ),
                        "Selection data",
                        html.Hr(),
                        html.Div(id='default-molecule3d-output')
                    ]), width=4),
        dbc.Col(mol_viewer_colorbar, width=4),
        dbc.Col(dbc.Card([dbc.CardBody(text_abreviation)
                    ]), width=4
            )

                ])
            ])



# app.layout = \
#     dbc.Row([dbc.Col(html.Div(html.B('text'), style={'height': '190px'},
#                 className='bg-danger'), width=4),
#              dbc.Col(html.Div([
#                 dbc.Row(dbc.Col(html.Div(html.B('overview'), className='bg-success')),
#                 dbc.Row([dbc.Col(html.Div(html.B('hist'), className='bg-success')),
#                          dbc.Col(html.Div(html.B('3dscatter'), className='bg-success'))]),
#                 dbc.Row([dbc.Col(html.Div([html.B('3D prot'), className='bg-success'])),
#                          dbc.Col(html.Div(html.B('rna'), className='bg-success')),
#                          dbc.Col(html.Div(html.B('text'), className='bg-success'))])
#                         )], width=8))])



def histogram(x_axis):
    fig = px.histogram(df, x=x_axis, color="clinvar_simple", marginal="rug")
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    fig.update_layout(
        height=250,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, visible=True, linecolor="gray", linewidth=5),
    )
    return fig

# Callback allows components to interact--------------------------------------------------------------------------------------------------
@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    Output(overview_title, 'children'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value'),
)
def update_overview_graph(column_name, y_axis_nucleotide):  # function arguments come from the component property of the Input
    print(column_name)
    print(y_axis_nucleotide)
    print(type(column_name))
    size_marker, height_grph, marker_symb = 8, 340, "square"
    limit = (0.5, 1.25)
    y_axis = [limit[0] + (limit[1]- limit[0])/2 for x in df['Genomic position']]
    yaxis_dict = dict(showgrid=False, visible=False)

    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        marker_symb, size_marker, y_axis, height_grph = "square", 8, 'alt_pos', 340
        yaxis_dict = dict(showgrid=False, visible=True, title='Nucleotide', tickvals=[0.5, 0.75, 1, 1.25],
                          ticktext=['T', 'G', 'C', 'A'])
    fig = px.scatter(data_frame=df,
                     x=df['Genomic position'],
                     y=y_axis,
                     height=height_grph,
                     # size="1/AC",
                     color=column_name,
                     custom_data=["SGE function classification",
                     "UKB function classification", 'Variant consequence', 'Clinvar annotation',
                                  'cohort_allele_count', 'var_name', 'Exon'],
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
        ref_plot = px.scatter(data_frame=df,
                              x=df['Genomic position'],
                              y=df['ref_pos'],
                              height=height_grph,
                              custom_data=["Genomic position", "ref", "Exon"])
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
                      line=dict(color='rgba(0,0,255, 0.5)', width=2),
                      fillcolor='rgba(0,0,255,0.15)', layer='below')
        if i + 1 in [23,10]:
            start = start +250
        if i + 1 in [ 17]:
            start = start +150
        fig.add_annotation(x=start, y=limit[1]+0.1,
                           text=texex,
                           showarrow=False,
                           font=dict(color='blue', size=8), textangle=270, xanchor=pos_xanchor, yanchor='bottom', bgcolor="white", opacity=1)

    # add intron
    for i in range(len(intron_list)):
        start, end = intron_list[i][0] - 0.5, intron_list[i][1] + 0.5
        fig.add_shape(type='line', x0=start, x1=end,
                      y0=limit[0] + (limit[1] - limit[0]) / 2, y1=limit[0] + (limit[1] - limit[0]) / 2,
                      line=dict(color='rgba(0,0,255, 0.5)', width=2), layer='below')

    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=5),
        yaxis=yaxis_dict,
        legend=dict(orientation='h', yanchor='top', y=1.5, xanchor='left', x=0, title="Annotation"))

    return fig, '### BRCA1 gene annotated with ' + column_name.lower()  # returned objects are assigned to the component property of the Output


# define a white to red gradient
colors = ['#FFFFFF', '#FF0000']
cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', colors)
@app.callback(
    Output('default-molecule3d-output', 'children'),
    Output(mol_viewer_colorbar, 'figure'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    # Get a color bar
    fig=px.scatter(df, x='rna.score.1', y='rna.score.2', color='minmax_neg_func_score', color_continuous_scale=['#FFFFFF', '#FF0000'],
               title="                  Comparison of RNA replicates")

    fig.update_traces(opacity=1)
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        yaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        height=450, width=550,
        title_font_size= 12,
        title_x=0.95,
        coloraxis_colorbar=dict(
            title="SGE fct score\n\n",
            # thicknessmode="pixels", thickness=50,
             lenmode="pixels", len=390,
            yanchor="top", y=1.3,x=-1.5,
            tickvals=[0,1],
            tickmode='array',
            ticktext=["Neutral", "LoF"],
            # ticks="outside", ticksuffix=" bills",
            # dtick=5
        ))

    if atom_ids is None or len(atom_ids) == 0:
        phr1='No amino acid has been selected. Click somewhere on the protein \
        structure to select an amino acid.'

        return 'No amino acid selected', fig
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


    if slct_data is None or slct_data == {'points': []}:
        fig2 = px.scatter_3d(df, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                            color='Exon', custom_data=["var_name",'Exon',
                                  'cohort_allele_count', 'cadd_score', 'minmax_neg_func_score'], title="Comparison between UKB allele count, CADD and SGE function score for all variants in BRCA1")
        fig2.update_traces(hovertemplate = "<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        fig2.update_layout(scene=dict(
            xaxis_title='CADD score',
            yaxis_title='SGE fct score',
            zaxis_title='1/UKB_ac',
            bgcolor='rgba(0,0,0,0)',
            xaxis_backgroundcolor= 'rgba(0,0,0,0)',
            yaxis_backgroundcolor='rgba(0,0,0,0)',
            zaxis_backgroundcolor='rgba(0,0,0,0)',
            xaxis_color='blue',
            xaxis_gridcolor='blue',
            yaxis_gridcolor='blue',
            zaxis_gridcolor='blue',
            xaxis=dict(showgrid=True),
            yaxis=dict(showgrid=True),
            zaxis=dict(showgrid=True)),
            legend=dict(orientation='v', yanchor='top', y=1.05, xanchor='left', x=0, title="Region"),
        height=800)

        return fig2, fig3, fig4, fig5

    if slct_data['points'] == []:
        fig2 = px.scatter_3d(title="Please select at least one variant")
        return fig2, fig3, fig4, fig5
    else:
        # print(hov_data['points'][0]['customdata'][0])
        print(f'selected data: {slct_data}')

        exons = [slct_data['points'][i]['customdata'][-1] for i in range(len(slct_data['points']))]
        var = [slct_data['points'][i]['customdata'][-2] for i in range(len(slct_data['points']))]

        print(set(exons))
        dff2 = df[df.var_name.isin(var)]
        fig2 = px.scatter_3d(dff2, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color='Exon',
                             custom_data=["var_name", 'Exon',
                                          'cohort_allele_count', 'cadd_score', 'minmax_neg_func_score'],
                             title="Comparison between UKB allele count, CADD and SGE function score for variant in "+str(set(exons)).replace("{", '').replace("'", '').replace("}", ''))
        fig2.update_traces(hovertemplate = "<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        fig2.update_layout(scene=dict(
            xaxis_title='CADD score',
            yaxis_title='SGE fct score',
            zaxis_title='1/UKB_ac'))

        return fig2, fig3, fig4, fig5





# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8054)
