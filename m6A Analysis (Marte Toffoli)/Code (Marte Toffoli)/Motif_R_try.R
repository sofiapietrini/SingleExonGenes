library("BiocManager")
library("Biostrings")
library(tidyr)
library(dplyr)
library(stringr)


#procede with analyzing the number of Motif that are modified on the total

ICs <- read.csv("/home/marte/M6A/M6A/ICs_with_length_mod.tsv", sep = "\t")
IGs <- read.csv("/home/marte/M6A/M6A/IGs_with_length_mod.tsv", sep = "\t")

stableIGs <- IGs$Gene.stable.ID.version
write.csv(stableIGs, "IGs_cod", quote = FALSE, row.names = FALSE)

stableICs <- ICs$Gene.stable.ID.version
write.csv(stableICs, "ICs_cod", quote = FALSE, row.names = FALSE)

# Read the sequences from the FASTA file
s_ig = readDNAStringSet("/home/marte/M6A/M6A/mart_exportIGs.txt")
s_ic = readDNAStringSet("/home/marte/M6A/M6A/mart_exportICs.txt")

#to check if the number of elements is the same, if not find a solution
#from DNAseqSet to dataframe
dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
s_ig_df <- as.data.frame(dss2df(s_ig)) #it works

#divide a column to have a new column containing just gene names
s_ig_df <- s_ig_df %>%
  separate_wider_delim(names, delim = "|", names = c("Gene Stable ID", "Transcript Stable ID", "Names","Chr"), names_sep = "\t")

#check if the dimension is ok
s_ig_df <- s_ig_df[(s_ig_df$`names	Gene Stable ID` %in% IGs$Gene.stable.ID),] #in s_ig_df 2 genes less
#IGs <- IGs[(IGs$GeneName %in% s_ig_df$`names	Names` & IGs$Chromosome.scaffold.name %in% s_ig_df$`names	Chr`),] #solve the issue. 2 genes are expendable



#do the same here
dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
s_ic_df <- as.data.frame(dss2df(s_ic)) #it works

s_ic_df <- s_ic_df %>%
  separate_wider_delim(names, delim = "|", names = c("Gene Stable ID", "Transcript Stable ID", "Names","Chr"), names_sep = "\t")

s_ic_df <- s_ic_df[(s_ic_df$`names	Gene Stable ID` %in% ICs$Gene.stable.ID),]



# Define your subsequence (motif) that you want to search for
DRACH <- c("AAACA", "AAACT", "AAACC","AGACA", "AGACC", "AGACT", "GAACA", "GAACT", "GAACC", "GGACA", "GGACC", "GGACT", "TAACA", "TAACC", "TAACT", "TGACA", "TGACC", "TGACT")
occurrencesIG <- matrix(NA, nrow = nrow(s_ig_df), ncol = length(DRACH))
occurrencesIC <- matrix(NA, nrow = nrow(s_ic_df), ncol = length(DRACH))


for (i in 1:length(DRACH)){
  # Convert subsequence to a DNAString object
  subseq <- DRACH[i]
  #subset <- as.data.frame(subseq)
  # Function to count subsequence occurrences in each sequence
  count_occurrences <- function(seq, subseq) {
    # Count how many times subseq appears in seq
    return(countPattern(subseq, seq))
  }

  # Apply this function to each sequence in the FASTA file
  occurrencesIG[,i] <- unname(sapply(s_ig_df$seq, count_occurrences, subseq))
  
  #occurrencesIC[,i] <- unname(sapply(s_ic_df$seq, count_occurrences, subseq))
  
}

write.csv(occurrencesIG, "/home/marte/M6A/IGs_o")
write.csv(occurrencesIC, "/home/marte/M6A/ICs_o")


IGs_o <- read.csv("/home/marte/M6A/IGs_o", sep = ",")
ICs_o <- read.csv("/home/marte/M6A/ICs_o", sep = ",")
sum_ig_o <- rowSums(IGs_o)
#sum_o <- t(t(sum_o))

RatioIG <- sum_ig_o/IGs$count

sum_ic_o <- rowSums(ICs_o)
#sum_o <- t(t(sum_o))
RatioIC <- sum_ic_o/ICs$count

length(RatioIG) <- length(RatioIC)
data <- data.frame( "multi-exon genes" = RatioIC, 
                    "single-exon genes" = RatioIG 
                    ) 
boxplot(data, na.action = NULL,
        ylab = "num. of mod. per num. of motifs", outline = FALSE, col = c("steelblue", "coral"))

# would calculate the "number_modifications/gene length" ratio, which seems directly proportional to the m6A density. 
# This could be expressed, for example, as the number of m6A modifications per 100 nucleotides, or per 1000 nucleotides.

IGs_ratio = (IGs$count * 100)/IGs$Transcript.length..including.UTRs.and.CDS.
ICs_ratio = (ICs$count * 100)/ICs$Transcript.length..including.UTRs.and.CDS.

length(IGs_ratio) <- length(ICs_ratio)

data <- data.frame( multi_exon_genes = ICs_ratio, 
                    single_exon_genes = IGs_ratio) 
boxplot(data, na.action = NULL,
        ylab = "num. of mod per 100 bps", outline = FALSE, notch = TRUE, col = c("steelblue", "coral"))
#add colors

wilcox.test(ICs_ratio, IGs_ratio) # since p-value is smaller than 0.05 we have to reject the H0, p1 != p2

summary(ICs_ratio)
summary(IGs_ratio) #mean is larger than the ICs' mean


write.csv(ICs$GeneName, "ICS_names", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.csv(IGs$GeneName, "ensemble_ig", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.csv(ICs$GeneName, "ensemble_ic", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

chrom <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","22","mt","x","y")
intersect(IGs$Chromosome.scaffold.name, chrom)
length(match(chrom, IGs$Chromosome.scaffold.name))

risultato <- sapply(chrom, function(x) sum(IGs$Chromosome.scaffold.name == x))
risultato <- sum(risultato)
#Take the differential expression with multi-exon genes and see if the multi-exon ones are more correlated to diseases
#than the single exon ones (in proportion)
#look if the number of M6A could have an influence on the result