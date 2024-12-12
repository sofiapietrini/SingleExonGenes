import pandas as pd

f = open("hsap_SEG.fa", "r")
lines = f.readlines()
lines = [i.strip('\n') for i in lines]
#print(lines)

list = []
name_tot = ""
role = ""
for i in lines:
    name_tot = ""
    role = ""
    if i.startswith(">"):
        i = i.split(" ")
        for j in range(0, 3):
            i.remove(i[0])
        #list.append(i)

        l_i = len(i)

        #to get the role name only
        for k in range(0, l_i - 5):
            role = role + f" {i[k]}"

        #print(role)
        name_tot = role
        for k in range(l_i - 5,l_i):
            #name_tot = f"{role}\t{i[5]}"
            name_tot = name_tot + f"\t{i[k]}"
        list.append(name_tot)
#print(list)

#table creation
with open("exon_file.tsv", "w") as out:
    out.write("function\tcod\tgene_name\ttype\tchr\torganism\n")
    for n in list:
        out.write(f"{n}\n")

values = pd.read_csv('exon_file.tsv', delimiter='\t')

#count how many times a gene is present for each function
gene = [riga[1] for riga in values]
gene = pd.DataFrame(values)

freq = gene['function'].value_counts()
        
df = pd.DataFrame(freq)
df.to_csv("output.tsv", sep='\t', index=True)

print(freq)


#data frame     
#2285
#print(list)