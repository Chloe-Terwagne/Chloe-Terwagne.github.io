"""
    Project :  variants intersection df for BRCA1
    Description: merge raw BRCA1 SGE data from G. Findlay et al., 2018 and the UKB allele table
    Name : Chloé Terwagne
    date : 19th Jan 2023
    Python version 3.10
"""
# IMPORT ---------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import datetime
from liftover import get_lifter
import matplotlib as plt


def clean_allele_table_ukb(df, gene):
    dict_rename_columnn = {'Consequence (Most severe by gene)': 'consequence', 'Variant ID': 'variant_id', 'Chr': 'chr',
                           'Position': 'position', 'RSID': 'rsid', 'Reference': 'ref_long', 'Alternate': 'alt_long',
                           'Type': 'type',
                           'Cohort AF': 'cohort_af', 'Cohort Allele Count': 'cohort_allele_count',
                           'GnomAD AF': 'gnomad_af'}

    dict_consequence_rename = {'["missense_variant:' + gene + '"]': 'missense',
                               '["stop_gained:' + gene + '"]': 'stop_gained',
                               '["start_lost:' + gene + '"]': 'start_lost',
                               '["splice_acceptor_variant:' + gene + '"]': 'splice_acceptor',
                               '["splice_donor_variant:' + gene + '"]': 'splice_donor',
                               '["synonymous_variant:' + gene + '"]': 'synonymous',
                               '["intron_variant:' + gene + '"]': 'intron',
                               '["splice_region_variant:' + gene + '"]': 'splice_region',
                               '["3_prime_UTR_variant:' + gene + '"]': "3_utr",
                               '["5_prime_UTR_variant:' + gene + '"]': "5_utr"}
    dict_assignation = {'missense': 'neutral',
                        'start_lost': 'lof', 'stop_gained': 'lof', 'splice_acceptor': 'lof', 'splice_donor': 'lof',
                        'synonymous': 'neutral', 'intron': 'neutral', 'splice_region': 'neutral',
                        "3_utr": 'neutral', "5_utr": 'neutral'}

    df.rename(dict_rename_columnn, axis=1, inplace=True)
    # for some reason the position is sometime incorrect i.e: 17_43124029_T_G as a position of 43124027 however
    # according to the rsid the true position is 43124029 reinforme position based on variand id
    new = df['variant_id'].str.strip()
    new = new.str.split('_', expand=True)
    df["position"] = new[1].astype('int64')
    # remove != SNP
    df = df[df.type == "SNP"]
    # clean ref and alt column
    df['alt'] = df['variant_id'].str.strip().str[-1]
    df['ref'] = df['variant_id'].str.strip().str[-3]
    df['consequence'].replace(dict_consequence_rename, inplace=True)
    # consequence assignation
    df["func_class"] = df['consequence'].map(dict_assignation)
    print(df.head())
    df = df[["variant_id", "chr", "position", "ref", "alt", "func_class", "consequence", "rsid", "cohort_af",
             "cohort_allele_count"]]
    return df


def liftover_df_hg19_to_hg38(df, col_chr_name, col_pos_name, inplace=True):
    """
    :param df: pandas dataframe to liftover
    :param col_chr_name: string, name of the hg19 chromosome column
    :param col_pos_name: string, name of the hg19 position column
    :param inplace: when True, chromosome and position column are replaced by the hg38 version. if False, new columns
    are created
    :return: the liftovered pandas dataframe
    """
    # Get position of all variants and initiate lists
    list_chr_hg19 = df[col_chr_name].to_list()
    list_pos_hg19 = df[col_pos_name].to_list()
    list_chr_hg38, list_pos_hg38 = [], []

    # liftover converter
    converter = get_lifter('hg19', 'hg38')
    for i in range(len(list_chr_hg19)):
        hg38 = converter[list_chr_hg19[i]][list_pos_hg19[i]]
        list_chr_hg38.append(int(hg38[0][0].replace('chr', '')))
        list_pos_hg38.append(int(hg38[0][1]))
    if inplace:
        df.rename({col_chr_name: 'chromosome_hg38', col_pos_name: 'position_hg38'}, axis=1, inplace=True)
    else:
        df.rename({col_chr_name: 'chromosome_hg19', col_pos_name: 'position_hg19'}, axis=1, inplace=True)
    df['chromosome_hg38'] = list_chr_hg38
    df['position_hg38'] = list_pos_hg38

    return df


# -- Displaying function


def print_msg_box(msg, indent=1, width=None, title=None):
    """Print message-box with optional title."""
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print(box)


def print_df_info(df, df_name, list_value_counts=[]):
    """Print basic dataframe information"""
    print_msg_box(df_name)
    print(df.head(3), '\n')
    print(df.shape, '\n')
    print("position variant interval: ", df["position"].min(), '-', df["position"].max(), '\n')
    if len(list_value_counts) > 0:
        for i in range(len(list_value_counts)):
            print(df[list_value_counts[i]].value_counts(), '\n')
    if len(list_value_counts) > 1:
        print(df.groupby(list_value_counts[0])[list_value_counts[1]].value_counts(), '\n')


# -- Plot function

def plot_position_var(position_list, name_df_list, color_list, alt_list, ref_list, list_position_conflict):
    """
    Warning: works only inside a same chromosome, plot variants from different source (i.e observed in UKB, tested in
    SGE or/and intersection)
    :param position_list: list of lists of the position of the variant, one list from each source
    :param name_df_list: list of lists of string, names of the source of variant
    :param color_list: list of lists of string, color to label each variant
    :param alt_list: list of lists of type of nucleotide variants (ATCG)
    :param alt_list: list of lists of type of nucleotide reference (ATCG)
    :param list_position_conflict:list[0] is a list of position where annotation conflit happen. list[1] correspond to
     the color desired for this conflit and list [2] bool that don't display the vertical conflict line when set to
     false.
    :return:
    """
    # Create a ATGC dict to spread on vertical axis
    ACGT_dict = {'A': 0.08, 'C': 0.025, 'G': -0.025, 'T': -0.08}

    # Figure size (width, height)
    plt.figure(figsize=(12, 6))
    y_axis = [i / 3 for i in range(len(position_list))]  # generate anchor point for the sources on the y axis
    plt.yticks(y_axis, name_df_list, rotation=90, va="center")  # adding label at anchor point position

    # plot ref and alt allele for each sources
    for i in range(len(name_df_list)):
        x = [j for j in position_list[i]]  # x genomic position
        y = [n + y_axis[i] for n in pd.Series(alt_list[i]).map(
            ACGT_dict)]  # y axis genomic position for alt (slightly shifted from the anchor point depending on the ACTG dict )
        y_ref = [n + y_axis[i] for n in pd.Series(ref_list[i]).map(
            ACGT_dict)]  # y axis genomic position for ref (slightly shifted from the anchor point depending on the ACTG dict )
        plt.scatter(x, y_ref, marker='s', alpha=1, edgecolors='black', color='white', zorder=2)
        plt.scatter(x, y, label=name_df_list[i], marker="s", alpha=1, zorder=3, color=color_list[i])

    # adjusting limits x-axis (for the horizontal line)
    min_pos, max_pos = min([min(i) for i in position_list]), max([max(i) for i in position_list])
    interval = max_pos - min_pos
    space = interval * 0.05
    plt.xlim(min_pos - space, max_pos + space)

    # Plotting gray horizontal lines and ATCG letter
    acgt_shit = [val for key, val in ACGT_dict.items()]
    acgt = [key for key, val in ACGT_dict.items()]
    for major in y_axis:
        for i in range(len(acgt_shit)):
            plt.hlines(major + acgt_shit[i], linewidth=1, color='grey', alpha=0.3, zorder=1,
                       xmin=min_pos - space, xmax=max_pos + space)
            plt.text(min_pos - space / 2, major + acgt_shit[i] - 0.005, acgt[i], fontsize=10)

    # adjusting limit y axis for the vertical line
    y_min, y_max = min(y_axis) - 0.09, max(y_axis) + 0.09
    if list_position_conflict[2]:
        for i in range(len(list_position_conflict[0])):
            plt.vlines(list(list_position_conflict[0])[i], linewidth=1.2, color=str(list(list_position_conflict[1])[i]),
                       alpha=1, zorder=1,
                       ymin=y_min, ymax=y_max)

    plt.xlabel('Genomic position')
    plt.ylabel('variants')
    plt.title("Variant positions")
    plt.tight_layout()
    # uncomment to save the fig
    # plt.savefig('variant_position_vhl_sge_intersection' + DATE + '.png', dpi=1200)
    plt.show()



pd.set_option('display.width', 600)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)

# CONSTANT -------------------------------------------------------------------------------------------------------------
DATE = str(datetime.date.today()).replace('-', '_')
GENE = "brca1"
PATH_AT = "/Users/terwagc/variant_prioritisation_data/input_df/allele_table_brca1_chr17_43044295_43125364.csv"
PATH_ASSAY = "/Users/terwagc/variant_prioritisation_data/input_df/subset_brca1_region/sge_paper_brac1_2018/Supplementary_Table_1_revised.csv"
FILE_OUTPUT = "/Users/terwagc/PycharmProjects/dataviz_brca1/df/merged_"+ GENE +"_sge_ukb_"+DATE+".csv"
VERBOSE = True

# MAIN -----------------------------------------------------------------------------------------------------------------

ukb_df = pd.read_csv(PATH_AT)
print("--ukb dataset loaded")
print(ukb_df.head(3))

sge_df = pd.read_csv(PATH_ASSAY)
print("--sge dataset loaded")
print(sge_df.head(143))

# SGE df ---------------------------------------------------------------------------------

sge_df = liftover_df_hg19_to_hg38(sge_df, "chromosome", "position (hg19)", inplace=True)
sge_df["variant_id"] = sge_df['chromosome_hg38'].astype(str) + '_' + sge_df['position_hg38'].astype(str) + '_' + \
                       sge_df['reference'] + '_' + sge_df['alt']
sge_df.rename({'chromosome_hg38': 'chr', 'reference': 'ref', 'function.score.mean': 'func_score',
               'phyloP (mammalian)': 'phylop_mam', 'position_hg38':'position',
               'func.class': 'func_class', 'CADD.score': 'cadd_score', 'gnomAD_AF': 'gnomad_af'},
              axis=1, inplace=True)
# rename consequence
name_consequence_sge = {'Missense': 'missense', 'Nonsense': 'nonsense',
                        'Canonical splice': 'canonical_splice',
                        'Synonymous': 'synonymous', 'Intronic': 'intron',
                        'Splice region': 'splice_region',
                        "5' UTR": "5_utr"}
sge_df['consequence'].replace(name_consequence_sge, inplace=True)
assignation_dict_sge = {'FUNC': 'neutral', 'LOF': 'lof', 'INT': 'intermediate'}
sge_df["func_class"] = sge_df['func_class'].map(assignation_dict_sge)
sge_df = sge_df[["variant_id", 'chr', 'position', 'ref', 'alt', 'func_score', 'func_class', 'consequence',
                 'cadd_score', 'phylop_mam', 'polyphen2', 'sift', 'clinvar', 'gnomad_af', "bravo_AF","aa_pos", 'aa_ref', 'aa_alt', 'protein_variant', "rna.score.1",  "rna.score.2"]]
print("-- sge df reformated")

# UKB df ---------------------------------------------------------------------------------
ukb_df = clean_allele_table_ukb(ukb_df,"BRCA1")
print("-- ukb allele table reformated")

if VERBOSE:
    print_df_info(sge_df, "SGE - BRCA1", ["func_class", "consequence"])
    print_df_info(ukb_df, "UKB - BRCA1", ["func_class", "consequence"])

# intersection UKB and sge by variant ID -------------------------------------------------------------------------------
intersec_sge_ukb_df = pd.merge(sge_df, ukb_df, on="variant_id", suffixes=('_sge', '_ukb'))

# Dropping identical columns
for col in ["ref", "alt", "chr", "position"]:
    if intersec_sge_ukb_df[col+'_sge'].astype(str).equals(intersec_sge_ukb_df[col+'_ukb'].astype(str)):
        intersec_sge_ukb_df[col] = intersec_sge_ukb_df[col+'_sge']
        intersec_sge_ukb_df.drop(columns=[col+'_ukb', col+'_sge'], inplace=True)

# is the function class in SGE different from the function class assigned by ukb?
intersec_sge_ukb_df["consequence_conflict"] = np.where(
    intersec_sge_ukb_df['consequence_ukb'] != intersec_sge_ukb_df['consequence_sge'], "yes", "no")

# is the function class in SGE different from the function class assigned by ukb?
intersec_sge_ukb_df["func_conflict"] = np.where(
    intersec_sge_ukb_df['func_class_ukb'] != intersec_sge_ukb_df['func_class_sge'],
    "yes", "no")
# How many conflit function have also conflict consequence?
intersec_sge_ukb_df["double_conflict"] = np.where((intersec_sge_ukb_df['func_conflict'] == "yes") & (
        intersec_sge_ukb_df['consequence_conflict'] == "yes"),
                                                  "yes", "no")
# add gray zone
intersec_sge_ukb_df["func_conflict"] = np.where(
    (intersec_sge_ukb_df['func_conflict'] == 'yes') & (intersec_sge_ukb_df['func_class_sge'] == 'intermediate'),
    "gray_zone", intersec_sge_ukb_df['func_conflict'])

print(intersec_sge_ukb_df.head())
# reordering column
intersec_sge_ukb_df = intersec_sge_ukb_df[
    ['variant_id', 'chr', 'position', 'ref', 'alt', 'func_score', 'func_class_sge', "func_class_ukb", "consequence_sge",
     "consequence_ukb", "func_conflict", "consequence_conflict", 'double_conflict', 'clinvar', 'cohort_allele_count',
     'cohort_af', 'rsid', 'cadd_score', 'phylop_mam', 'polyphen2','sift', "bravo_AF","aa_pos", 'aa_ref', 'aa_alt', 'protein_variant', "rna.score.1",  "rna.score.2"]]


print("-- dataset merged ----------------------------------------------------------------")
if VERBOSE:
    print_msg_box("Merged dataset")
    print(intersec_sge_ukb_df.head(5))
    print("--sge-ukb dataset nb variant: ", intersec_sge_ukb_df.shape[0])
    print("--sge-ukb dataset nb variant unique position:", len(intersec_sge_ukb_df["position"].unique()))
    print("--sge-ukb dataset position interval variant: ", intersec_sge_ukb_df["position"].min(), '-',
          intersec_sge_ukb_df["position"].max())

    print("\n--sge-ukb dataset consequence_sge classes counts:\n", intersec_sge_ukb_df['consequence_sge'].value_counts())
    print("\n--sge-ukb dataset consequence_ukb classes counts:\n", intersec_sge_ukb_df['consequence_ukb'].value_counts())
    print("\n--Do sge-ukb dataset has conflict variant consequence attributions?:\n",
          intersec_sge_ukb_df['consequence_conflict'].value_counts())

    print("\n\n--sge-ukb dataset functional classes counts:\n", intersec_sge_ukb_df['func_class_sge'].value_counts())
    print("\n--sge-ukb dataset functional classes counts:\n",
          intersec_sge_ukb_df['func_class_ukb'].value_counts())
    print("\n--Do sge-ukb dataset has conflict variant functional score/annotation attributions?:\n",
          intersec_sge_ukb_df['func_conflict'].value_counts())
intersec_sge_ukb_df.to_csv(FILE_OUTPUT, index=False)
print(intersec_sge_ukb_df.describe())
print("-- saved as", FILE_OUTPUT)
