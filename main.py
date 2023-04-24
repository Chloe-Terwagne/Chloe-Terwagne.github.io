# If you prefer to run the code online instead of on your computer click:
# https://github.com/Coding-with-Adam/Dash-by-Plotly#execute-code-in-browser
import numpy as np
from dash import Dash, dcc, Output, Input  # pip install dash
import dash_bootstrap_components as dbc  # pip install dash-bootstrap-components
import plotly.express as px
import pandas as pd  # pip install pandas
from dash import dcc, html, Output, Input
from protein_folding import create_style
import dash_bio as dashbio
from dash import html
from dash_bio.utils import PdbParser

pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)


def adding_cols(df, exons):
    print(df.head(50))
    df = df.sort_values("position")
    df = df.rename(columns={"position": "genomic position"})
    df['var_index'] = [x for x in range(len(df['genomic position']))]
    # normalize function score
    df['minmax_neg_func_score'] = [-a for a in df['func_score'].to_list()]
    df['minmax_neg_func_score'] = (df['minmax_neg_func_score'] - df[
        'minmax_neg_func_score'].min()) / (df['minmax_neg_func_score'].max() - df[
        'minmax_neg_func_score'].min())
    df['clinvar_simple'] = "absent / no clear interpretation"
    df['clinvar_simple']=df['clinvar'].map({"Likely pathogenic":"Pathogenic", 'Pathogenic':'Pathogenic',"Pathogenic/Likely pathogenic":"Pathogenic", "Likely benign":"Benign","Benign":'Benign'})

    # Get size depend on AC
    df['1/AC'] = [1 / x for x in df['cohort_allele_count'].to_list()]
    # Get Y axis based alt nucleotide var
    df['alt_pos'] = [x for x in df['alt'].map({"A": 2, "C": 1.5, "G": 1, "T": 0.5})]
    df['ref_pos'] = [x for x in df['ref'].map({"A": 2, "C": 1.5, "G": 1, "T": 0.5})]
    df['var_name'] = "chr" + df['chr'].astype(str) + " " + df['genomic position'].astype(str) + " " + df['ref'].astype(
        str) + ">" + df['alt'].astype(str)
    df["exon"]= "intron"
    for i in range(len(exons)):
        df["exon"] = np.where((df["genomic position"]> exons[i][1]) & (df["genomic position"]<exons[i][0]), "exon "+str(i+1), df["exon"])

    print(df.head())
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
#3D parser
parser = PdbParser('/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/AF-P38398-F1-model_v4.pdb')
# from https://alphafold.ebi.ac.uk/entry/P38398
data = parser.mol3d_data()

styles = create_style(
    df, data['atoms'], visualization_type='cartoon', color_element='residue_score'
) # TODO add the column of score we want to add
print("---------------")

overview_title = dcc.Markdown(children='')
overview_display = dcc.RadioItems(options=["collapsed\t", "expand nucleotide type"], value='collapsed\t')
overview_dropdown = dcc.Dropdown(options=['clinvar', 'consequence_ukb', "func_class_sge", "func_class_ukb"],
                                 value='consequence_ukb',  # initial value displayed when page first loads
                                 clearable=False)
overview_graph = dcc.Graph(figure={}, config={
                      'staticPlot': False,     # True, False
                      'scrollZoom': False,      # True, False
                      'doubleClick': 'reset',  # 'reset', 'autosize' or 'reset+autosize', False
                      'showTips': True,       # True, False
                      'displayModeBar': True,  # True, False, 'hover'
                      'watermark': False,
                      # 'modeBarButtonsToRemove': ['pan2d','select2d'],
                        }, selectedData=None)

# select Graph
three_d_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': True,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})
three_d_title = dcc.Markdown(children='all variant')
#
clinvar_hist_graph_sge = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})
clinvar_hist_graph_ac = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})
clinvar_hist_graph_cadd = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})

# rna_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})
# shankley_d_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False,'doubleClick': 'reset','showTips': True,'displayModeBar': True,'watermark': False})

# Customize your own Layout--------------------------------------------------------------------------------------------------
app.layout = html.Div([

    html.H1("Exploration of BRCA1", style={'text-align': 'center'}),

    overview_title,
    overview_display,
    overview_dropdown,
    overview_graph,
    html.Br(),
    three_d_graph,

    html.Div([
        dashbio.Molecule3dViewer(
            id='dashbio-default-molecule3d',
            modelData=data,
            styles=styles,
            selectionType='residue',
        ),
        "Selection data",
        html.Hr(),
        html.Div(id='default-molecule3d-output')
    ]),

    clinvar_hist_graph_sge,
    clinvar_hist_graph_ac,
    clinvar_hist_graph_cadd
])

def histogram(x_axis):
    fig = px.histogram(df, x=x_axis, color="clinvar_simple", marginal="rug")
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, visible=True, linecolor="gray", linewidth=5),
    )
    return fig

# AUTO GENERATED FILE - DO NOT EDIT




# Callback allows components to interact--------------------------------------------------------------------------------------------------
@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    Output(overview_title, 'children'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value')
)
def update_overview_graph(column_name,
                          y_axis_nucleotide):  # function arguments come from the component property of the Input
    print(column_name)
    print(y_axis_nucleotide)
    print(type(column_name))
    size_marker, height_grph, marker_symb = 12, 200, "line-ns-open"
    y_axis = [2 for x in df['genomic position']]
    limit = (1.5, 2.5)
    yaxis_dict = dict(showgrid=False, visible=False)

    print(limit[0] + (limit[1] - limit[0]) / 2)
    if y_axis_nucleotide == "expand nucleotide type":
        marker_symb, size_marker, y_axis, height_grph = "square", 8, 'alt_pos', 300
        limit = (0, 2.5)
        yaxis_dict = dict(showgrid=False, visible=True, title='Nucleotide', tickvals=[0.5, 1, 1.5, 2],
                          ticktext=['T', 'G', 'C', 'A'])
    fig = px.scatter(data_frame=df,
                     x=df['genomic position'],
                     y=y_axis,
                     height=height_grph,
                     # size="1/AC",
                     color=column_name,
                     custom_data=['func_class_sge', 'func_class_ukb', 'consequence_sge', 'clinvar',
                                  'cohort_allele_count', 'var_name', 'exon'],
                     )

    fig.update_traces(marker=dict(size=size_marker, symbol=marker_symb), selector=dict(mode="markers"),
                      hovertemplate="<br>".join([
                          "<b>%{customdata[5]}</b>",
                          #"Genomic position: %{x}",
                          "SGE classification: %{customdata[0]}",
                          "UKB classification: %{customdata[1]}",
                          "Clinvar classication: %{customdata[3]}",
                          "Consequence: %{customdata[2]}",
                          "UKB allele count: %{customdata[4]}",
                      ]))

    # add reference variant
    if y_axis_nucleotide == "expand nucleotide type":
        ref_plot = px.scatter(data_frame=df,
                              x=df['genomic position'],
                              y=df['ref_pos'],
                              height=height_grph,
                              hover_name=["Reference allele" for var in df["variant_id"]],
                              custom_data=['func_class_sge', 'func_class_ukb', 'consequence_sge', 'clinvar',
                                           'cohort_allele_count', 'var_name', 'exon'])
        ref_plot.update_traces(marker=dict(size=size_marker, symbol="square-open"), selector=dict(mode="markers"),
                               name='ref')
        # fig.update_yaxes(range=[limit[0]-1, limit[1]])
        fig.add_trace(ref_plot.data[0])

    # add exon
    for i in range(len(exon_list)):
        start, end = exon_list[i][1] - 0.5, exon_list[i][0] + 0.5
        fig.add_shape(type="rect",
                      x0=start, y0=limit[0], x1=end, y1=limit[1],
                      line=dict(color='rgba(0,0,255, 0.5)', width=2),
                      fillcolor='rgba(0,0,255,0.15)', layer='below')
        fig.add_annotation(x=start, y=limit[1],
                           text="exon" + str(i + 1),
                           showarrow=False,
                           font=dict(color='blue', size=8), textangle=290, xanchor='left', yanchor='bottom')
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
    )

    return fig, '## Overview BRCA1 variant ' + column_name + " annotated"  # returned objects are assigned to the component property of the Output


@app.callback(
    Output('default-molecule3d-output', 'children'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    if atom_ids is None or len(atom_ids) == 0:
        return 'No atom has been selected. Click somewhere on the molecular \
        structure to select an atom.'
    return [html.Div([ #TODO change format to show fct score
        #html.Div('Element: {}'.format(data['atoms'][atm]['elem'])),
        #html.Div('Chain: {}'.format(data['atoms'][atm]['chain'])),
        html.Div('Residue name: {}'.format(data['atoms'][atm]['residue_name'])),
        html.Div('Residue name: {}'.format(df.loc[df['aa_pos'] == data['atoms'][atm]['residue_index'], 'func_score'])),

        html.Br()
    ]) for atm in atom_ids]

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
                            color='exon', custom_data=['func_class_sge', 'func_class_ukb', 'consequence_sge', 'clinvar',
                                                       'cohort_allele_count', 'var_name'], title="Rareness vs function score vs computational score for all var")
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
                             color='exon',
                             custom_data=['func_class_sge', 'func_class_ukb', 'consequence_sge', 'clinvar',
                                          'cohort_allele_count', 'var_name'],
                             title="Rareness vs function score vs computational score for var in "+str(set(exons)).replace("{", '').replace("'", '').replace("}", ''))
        return fig2, fig3, fig4, fig5





# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8054)
