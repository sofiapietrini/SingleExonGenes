import pandas as pd

# Read the FASTA file
with open("hsap_SEG.fa", "r") as f:
    lines = f.readlines()
    lines = [i.strip('\n') for i in lines]

# List to store the formatted rows
formatted_rows = []

# Process each line
for i in lines:
    if i.startswith(">"):
        # Split the header on spaces
        i = i.split(" ")
        
        # Remove the first element (assuming it's not needed)
        i.pop(0)

        l_i = len(i)

        # Initialize a list to collect fields for the row
        row = []

        # Collect fields based on your needs
        if l_i > 10:  # Ensure there are enough elements
            rows = row.extend(i[0:4])
            rows = row.extend(i[l_i - 5 : l_i -1])
            #row1 = sorted(row)
            #print(row1)
        formatted_rows.append(row)
            

# Write to the TSV file
with open("SingleExon_genes.tsv", "w") as out:
    # Write the header
    out.write("Strand\tlocus1\tlocus2\tcode\tcode2\tGene name\tSEG\tChr\n")
    
    # Write each row
    for row in formatted_rows:
        out.write("\t".join(row) + "\n")
        #print(out)