
ICs <- read.csv("/home/marte/M6A/M6A/ICs_with_length_mod.tsv", sep = "\t")
IGs <- read.csv("/home/marte/M6A/M6A/IGs_with_length_mod.tsv", sep = "\t")

Cs_mean <- mean(ICs_n)
IGs_mean <- mean(IGs_n)

ICs_median <- median(ICs_n)
IGs_median <- median(IGs_n)

ICs_IQR <- IQR(ICs_n)
IGs_IQR <- IQR(IGs_n)

ICs_quantile <- quantile(ICs_n, prob=c(.25,.5,.75), type=1)
IGs_quantile <- quantile(IGs_n, prob=c(.25,.5,.75), type=1)


ICs_Length <- ICs$Length..Mod 
IGs_Length <- IGs$Length..Mod[1:(length(IGs$Length..Mod)-1)] 
#padding with NA
length(IGs_Length) <- length(ICs_Length)


data <- data.frame( ICs_Length, 
                    IGs_Length
                    ) 
boxplot(data)


