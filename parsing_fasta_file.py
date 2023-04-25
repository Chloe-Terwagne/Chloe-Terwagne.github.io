from Bio import SeqIO
import pandas as pd

# specify the path to your FASTA file
fasta_file = "/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/BRCA1_38.fna"

# create an empty dataframe with the desired column names
df = pd.DataFrame(columns=["chromosome", "position", "nucleotide"])

# parse the FASTA file and extract the relevant information
for record in SeqIO.parse(fasta_file, "fasta"):
    chromosome = record.id.split()[0]  # extract only the chromosome name from the header
    print(chromosome)
    sequence = str(record.seq)  # extract the nucleotide sequence
    for i, nucleotide in enumerate(sequence):
        position = 43170327 - i  # adjust for 0-indexing
        if i%1000==0:
            print(position)
        # append the information to the dataframe as a new row
        df_new_row = pd.DataFrame({"chromosome": [chromosome], "position": [position], "nucleotide": [nucleotide]})
        #print(df_new_row)
        df = pd.concat([df, df_new_row])
print(df.head())
df.to_csv( "/Users/terwagc/PycharmProjects/dataviz_brca1/Chloe-Terwagne.github.io/df/brca1_refseq.csv")
