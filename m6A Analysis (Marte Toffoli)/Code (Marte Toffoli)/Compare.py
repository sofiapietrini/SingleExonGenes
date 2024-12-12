import pandas as pd
from Bio import Align
from Bio.Seq import Seq

# Read the TSV files into DataFrames
df1 = pd.read_csv('seq.tsv', sep='\t')  
df2 = pd.read_csv('gene_sequences_with_ensembl_ids.csv', sep=',')
df3 = pd.read_csv('exon_file.tsv', sep='\t')

# Assuming you want to check if values in a specific column of df2 are in a specific column of df1
column_name_df2 = 'sequence'
column_name_df1 = 'sequence'

column2 = 'gene_name'
column3 = 'gene_name'
#if column_name_df1 in df1.columns and column_name_df2 in df2.columns:
# Find values in df2 that are present in df1
present_seq = df2[df2[column_name_df2].isin(df1[column_name_df1])]
present_names = df2[df2[column2].isin(df3[column3])]

present_seq.to_csv('present_seq_comp.tsv', sep='\t', index=False, header=True)
present_names.to_csv('present_names_comp.tsv', sep='\t', index=False, header=True)

def To_protein (dnastring):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
    protein = []
    dnastring = dnastring.replace("\n", "")
    end = len(dnastring) - (len(dnastring) %3) - 1
    for i in range(0,end,3):
        codon = dnastring[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("N")
    return "".join(protein)

# Read the TSV file

# Extract sequences from the DataFrame
sequences1 = df1['sequence'].tolist() 
sequences2 = df2['sequence'].tolist()

# Compare each sequence with every other sequence
for i in range(len(sequences1)):
    for j in range(len(sequences2)):  # Start from i + 1 to avoid duplicate comparisons
        seq1 = Seq(sequences1[i])
        #print(seq1)
        seq2 = To_protein(Seq(sequences2[j]))
        #print(seq2)

        aligner = Align.PairwiseAligner()

        # Set scoring parameters (optional)
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.gap_score = -2

        # Perform global alignment
        alignment = aligner.align(seq1, seq2)
        

        # Check if the sequences are identical
        if seq1 == seq2:
            #print(f"Row {i + 1} and Row {j + 1}: The sequences are identical.")
            # Print the best alignment
            print(f"Alignment between row {i + 1} and row {j + 1}:")
            for aln in alignment:
                print(aln)
        #else:
            #print(f"Row {i + 1} and Row {j + 1}: The sequences are not identical.")
