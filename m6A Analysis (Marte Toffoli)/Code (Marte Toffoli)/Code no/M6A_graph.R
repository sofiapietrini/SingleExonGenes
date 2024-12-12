
ICs <- read.csv("/home/marte/M6A/M6A/ICs_with_length_mod.tsv", sep = "\t")
IGs <- read.csv("/home/marte/M6A/M6A/IGs_with_length_mod.tsv", sep = "\t")

min_max_normalization <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

ICs_n = min_max_normalization(ICs$Length..Mod)
IGs_n = min_max_normalization(IGs$Length..Mod)

#barplot(df2$Length..Mod[1:20],names.arg = df2$GeneName[1:20], legend.text = "IGs")
#barplot(df1$Length..Mod[1:20],names.arg = df1$GeneName[1:20], legend.text = "ICs")

#IGs_table <- prop.table(table(df2$Length..Mod))
#ICs_table <- prop.table(table(df1$Length..Mod))

ICs_mean <- mean(ICs_n)
IGs_mean <- mean(IGs_n)

ICs_IQR <- IQR(ICs_n)
IGs_IQR <- IQR(IGs_n)

ICs_quantile <- quantile(ICs_n, prob=c(.25,.5,.75), type=1)
IGs_quantile <- quantile(IGs_n, prob=c(.25,.5,.75), type=1)

boxplot(ICs_n)

boxplot(IGs_n)

#procede with analyzing the number of Motif that are modified on the total