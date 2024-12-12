import pandas as pd
import scipy.stats as stats
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message=".*longdouble.*")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import gseapy as gp
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import mannwhitneyu
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial import ConvexHull
from matplotlib.patches import Ellipse
import gseapy as gp

turnover = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Bladder", 
    "Colon_Sigmoid", "Colon_Transverse", "Colon_Transverse_Mixed_Cell", 
    "Colon_Transverse_Mucosa", "Esophagus_Gastroesophageal_Junction", 
    "Esophagus_Mucosa", "Liver", "Liver_Hepatocyte", "Liver_Mixed_Cell", 
    "Liver_Portal_Tract", "Lung", "Skin_Not_Sun_Exposed_Suprapubic", 
    "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", 
    "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate", 
    "Small_Intestine_Terminal_Ileum_Mixed_Cell", "Spleen", "Stomach", 
    "Stomach_Mixed_Cell", "Stomach_Mucosa", "Testis", "Thyroid", "Whole_Blood"
]

non_turnover = [
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", 
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
    "Brain_Hippocampus", "Brain_Hypothalamus", 
    "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", 
    "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", 
    "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", 
    "Kidney_Medulla", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", 
    "Pancreas", "Pancreas_Acini", "Pancreas_Islets", "Pancreas_Mixed_Cell", 
    "Prostate", "Vagina",  "Colon_Transverse_Muscularis", "Stomach_Muscularis"
]

uncertain = [
    "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", 
    "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", 
    "Cells_EBV-transformed_lymphocytes", "Cervix_Ectocervix", 
    "Cervix_Endocervix", "Esophagus_Muscularis", "Fallopian_Tube", 
    "Minor_Salivary_Gland", "Pituitary", "Uterus"
]

df=pd.read_csv("AllGene1.tsv", sep="\t" )
remove=[]
for index, row in df.iterrows():
    somma = row[2:70].sum()
    if somma==0:
        remove.append(row.iloc[1])
    else:
        for col, value in row[2:70].items():
            df.at[index, col]=float(value)/float(somma)#necessario perchè alcuni hanno numeri assurdi
df = df[~df['Gene'].isin(remove)]
df = df.drop(columns=uncertain)
df['meanTurn'] = np.nan
df['meanNonTurn'] = np.nan
df['diff'] = np.nan
df["def"]=np.nan
df['def'] = df['def'].astype('object')
#df['def'] = df['def'].astype('object')
soglia=0.95
for index,row in df.iterrows():
    valT=[]
    valN=[]
    for col_name, value in row[2:].items():
        #print(col_name)
        if col_name in turnover:
            valT.append(value)
        elif col_name in non_turnover:
            valN.append(value)
    meanT=sum(valT)/len(valT)
    meanN=sum(valN)/len(valN)
    df.at[index, 'meanTurn']=meanT
    df.at[index, 'meanNonTurn']=meanN
    df.at[index, 'diff']= meanT-meanN
    if meanT>=soglia/len(valT):
        df.at[index, 'def']="Turn"
    elif meanN>=soglia/len(valN):
            df.at[index, 'def']="Nturn"
    else:
        df.at[index, 'def']="noSig"

df2=df[["diff","def"]]
mmax=max(df["diff"])
mmin=min(df["diff"])
print(len(df))
# Creazione dell'istogramma
sns.histplot(data=df2, x='diff', hue='def', kde=False, bins=24, multiple='stack', palette='Set2')
plt.xlim(mmin,mmax)
# Aggiungi titolo ed etichette
plt.xlabel('∑Turn/N°Turn - ∑NotTurn/N°NotTurn', fontsize=14)
plt.ylabel('Number of Gene', fontsize=14)

# Mostra il grafico
plt.savefig("hist.pdf", format="pdf")
plt.close()

#dff=df[(df["diff"]==1)]

df=df.sort_values(by='diff', ascending=True)
df["Gene"].to_csv("prova.tsv", sep="\t")
T1=df[df['def']=='Turn']
T=len(T1)
T1["Gene"].to_csv("TurnOver.tsv", sep="\t")
N=len(df[df['def']=='Nturn'])
N1=df[df['def']=='Nturn']
T1["Gene"].to_csv("TurnOver1.tsv", sep="\t")
N1["Gene"].to_csv("NTurnOver1.tsv", sep="\t")
tot=len(df)
print(f"Turn={T} NonTurn={N} %T={T/tot} %N={N/tot}")
