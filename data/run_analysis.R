library(ggplot2)
library(extrafont)
library(seqinr)
library(Biostrings)
library(biomaRt)
library(plyr)
library(seqinr)

install.packages("devtools")

#this is the newest version of Mario dos Reis package: https://github.com/mariodosreis/tai 
devtools::install_github("mariodosreis/tai")
require(tAI)



cbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_science <- function (base_size = 12, base_family = "Arial Rounded MT Bold") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

#read file human_transcripts.txt, contains Accession number, transcript length and 3'UTR length
human_transcripts <- read.table("human_transcripts.txt")
#read file for miR-155 transfection footprint data
miR_155_fp <- read.delim("mir155_32hr_fp.txt", stringsAsFactors=FALSE)
#read file for mock transfection footprint data
mock_fp <- read.delim("mock_32hr_fp.txt", stringsAsFactors=FALSE)
#merge miR-155 and mock fp data
miR_155_mock <- merge(miR_155_fp, mock_fp, by.x="Gene_Name", by.y="Gene_Name")
#make new column for log2 fold change in fp
miR_155_mock$fold <- log2(miR_155_mock$Expression_Level_.rpkM..x/miR_155_mock$Expression_Level_.rpkM..y)
#remove NAs, Inf and -Inf from fold
miR_155_mock <- na.omit(miR_155_mock)
miR_155_mock <- subset(miR_155_mock, fold >-Inf)
miR_155_mock <- subset(miR_155_mock, fold <Inf)
#read file for mock transfection RNAseq data
mock_mRNA <- read.delim("mRNA_mock_32hr.txt", stringsAsFactors=FALSE)
#merge miR-155/mock fp data with mock RNAseq data
miR_155_TE <- merge(miR_155_mock, mock_mRNA, by.x="Gene_Name", by.y="Gene_Name")
#calculate TE (fpsubmock/RNAseqsubmock) and make a new column, TE
miR_155_TE$TE <- miR_155_TE$Expression_Level_.rpkM..y/miR_155_TE$Expression_Level_.rpkM.
#remove Inf from TE
miR_155_TE <- subset(miR_155_TE, TE <Inf)
#read file for miR-155 target identities from Targetscan
miR_155_TS <- read.delim("mir-155_TS.txt")
#merge miR-155/mock/TE data with miR-155 target identities
miR_155_targets <- merge(miR_155_TE, miR_155_TS, by.x="Gene_Name", by.y="Ortholog.of.target.gene")
miR_155_targets <- subset(miR_155_targets, TE >0)

ggplot(miR_155_targets, aes(x=log2(TE), y= fold)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  ylab("RPF Fold Change (log2)") +
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  xlab("TE (log2)") +
  scale_y_continuous(limits = c(-4, 6), expand = c(0, 0)) +
  geom_text(family="Arial", fontface="italic", x = -2.5, y = 5, label = "r = -0.4953742, p = 2.2e-16") +
  theme_science()
ggsave("fold_logTE_lm.tiff",  width = 5, height = 3)
ggplot(miR_155_targets, aes(x=TE, y= fold)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  ylab("RPF Fold Change (log2)") +
  xlab("TE") +
  scale_y_continuous(limits = c(-4, 6), expand = c(0, 0)) +
  geom_text(family="Arial", fontface="italic", x = 7.5, y = 4, label = "r = -0.3224535, p = 9.382e-11") +
  theme_science()
ggsave("fold_TE_lm.tiff",  width = 5, height = 3)

#define x as the median of TE
x = median(miR_155_targets$TE)
#make new column to describe TE as High or Low, above or below median
miR_155_targets$TE.high.low <- ifelse(miR_155_targets$TE <x, "Low", "High")
#KS test for fold repression
ks.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="High"], miR_155_targets$fold[miR_155_targets$TE.high.low=="Low"])
ks.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="High"], miR_155_targets$fold)
ks.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="Low"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="High"], miR_155_targets$fold[miR_155_targets$TE.high.low=="Low"])
ansari.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="High"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE.high.low=="Low"], miR_155_targets$fold)


#make an ECDF plot of fold for all miR-155 targets and those with High or Low TE
ggplot(miR_155_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets, aes(x=fold, group=TE.high.low, colour=TE.high.low), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-155_32_TE.tiff", width = 5, height = 3)

median(miR_155_targets$fold[miR_155_targets$TE.high.low=="High"])
median(miR_155_targets$fold[miR_155_targets$TE.high.low=="Low"])
median(miR_155_targets$fold)

#define quartile values TE
q1_TE_32 <- quantile(miR_155_targets$TE, 0.25)
q2_TE_32 <- quantile(miR_155_targets$TE, 0.5)
q3_TE_32 <- quantile(miR_155_targets$TE, 0.75)
q4_TE_32 <- quantile(miR_155_targets$TE, 1)
#function to bin TE by quartile
TE <- function(x) { 
  if(x <=q1_TE_32) y <- "Low"
  if(x >=q1_TE_32 & x <=q2_TE_32) y <- "Med.Low"
  if(x >=q2_TE_32 & x <=q3_TE_32) y <- "Med.High"
  if(x >=q3_TE_32) y <- "High"
  return(y)
}
#applies above function to new column
miR_155_targets$TE_q <- sapply(miR_155_targets$TE,TE)
#make an ECDF plot of fold for all miR-155 targets and those binned by TE quartile
ggplot(miR_155_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets, aes(x=fold, group=TE_q, colour=TE_q), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Med.High", "Med.Low", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-155_32_TEq.tiff", width = 5, height = 3)

#KS test for fold repression, TE by quartiles compared to all
ks.test(miR_155_targets$fold[miR_155_targets$TE_q=="High"], miR_155_targets$fold)
ks.test(miR_155_targets$fold[miR_155_targets$TE_q=="Med.High"], miR_155_targets$fold)
ks.test(miR_155_targets$fold[miR_155_targets$TE_q=="Med.Low"], miR_155_targets$fold)
ks.test(miR_155_targets$fold[miR_155_targets$TE_q=="Low"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE_q=="High"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE_q=="Med.High"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE_q=="Med.Low"], miR_155_targets$fold)
ansari.test(miR_155_targets$fold[miR_155_targets$TE_q=="Low"], miR_155_targets$fold)


#find the median fold change for messages with each TE quartile
median(miR_155_targets$fold[miR_155_targets$TE_q=="High"])
median(miR_155_targets$fold[miR_155_targets$TE_q=="Med.High"])
median(miR_155_targets$fold[miR_155_targets$TE_q=="Med.Low"])
median(miR_155_targets$fold[miR_155_targets$TE_q=="Low"])
median(miR_155_targets$fold)

#make a new column that sums the number of conserved and poorly conserved sites
miR_155_targets$Number_of_Sites <- miR_155_targets$Conserved.sites.total + miR_155_targets$Poorly.conserved.sites.total
#make a box plot of fold change grouped by TE and binned by number of sites
ggplot(subset(miR_155_targets, Number_of_Sites <3), aes(y = fold, x = as.factor(Number_of_Sites), fill = TE.high.low)) + 
  geom_boxplot() +
  scale_fill_manual(name = "TE", values = c("#E69F00", "#56B4E9")) +
  ylab("RPF Fold Change (log2)") +
  xlab("Number of miR-155 Sites") +
  theme_science()
#save above plot
ggsave("miR-155_32_box.tiff", width = 5, height = 3)


#merge miR-155/mock/TE/targets with human_transcripts
miR_155_targets_UTR <- merge(miR_155_targets, human_transcripts, by.x="RefSeq_Accession.x", by.y="RefSeq.mRNA")
#deine the median of 3'UTR length and make a new column for 3'UTR length above and below median
c = median(miR_155_targets_UTR$UTR_length)
miR_155_targets_UTR$UTR_median <- ifelse(miR_155_targets_UTR$UTR_length <c, "Short", "Long")
#KS test for fold change difference between lon and short 3'UTR messages
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$UTR_median=="Long"], miR_155_targets_UTR$fold[miR_155_targets_UTR$UTR_median=="Short"])
#ECDF plot of 3'UTR length and fold change
ggplot(miR_155_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_UTR, aes(x=fold, group=UTR_median, colour=UTR_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "3'UTR Length", values=cbPalette, breaks=c("Long", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_UTR.tiff", width = 5, height = 3)
#define 3'UTR length quartiles
q1_32 <- quantile(miR_155_targets_UTR$UTR_length, 0.25)
q2_32 <- quantile(miR_155_targets_UTR$UTR_length, 0.5)
q3_32 <- quantile(miR_155_targets_UTR$UTR_length, 0.75)
q4_32 <- quantile(miR_155_targets_UTR$UTR_length, 1)
#function for 3'UTR length to name transcripts in each quartile
Length <- function(x) { 
  if(x <=q1_32) y <- "Short"
  if(x >=q1_32 & x <=q2_32) y <- "Med.Short"
  if(x >=q2_32 & x <=q3_32) y <- "Med.Long"
  if(x >=q3_32) y <- "Long"
  return(y)
}
#applies the above function
miR_155_targets_UTR$Length_q <- sapply(miR_155_targets_UTR$UTR_length, Length)
#ECDF plot of 3'UTR length and fold change
ggplot(miR_155_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_UTR, aes(x=fold, group=Length_q, colour=Length_q), stat="ecdf", size=1) +
  scale_colour_manual(name= "3'UTR Length", values=cbPalette, breaks=c("Long", "Med.Long", "Med.Short", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_UTRq.tiff", width = 5, height = 3)
#ks tests for fold change between 3'UTR length quartiles
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_q=="Long"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_q=="Med.Long"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_q=="Med.Short"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_q=="Short"], miR_155_targets_UTR$fold)

#deine the median of transcript length and make a new column for transcript length above and below median
d = median(miR_155_targets_UTR$Transcript_length)
miR_155_targets_UTR$Transcript_median <- ifelse(miR_155_targets_UTR$Transcript_length <d, "Short", "Long")
#KS test for fold change difference between long and short messages
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Transcript_median=="Long"], miR_155_targets_UTR$fold[miR_155_targets_UTR$Transcript_median=="Short"])
#ECDF plot of length and fold change
ggplot(miR_155_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_UTR, aes(x=fold, group=Transcript_median, colour=Transcript_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "Transcript Length", values=cbPalette, breaks=c("Long", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_length.tiff", width = 5, height = 3)
#define length quartiles
qL1 <- quantile(miR_155_targets_UTR$Transcript_length, 0.25)
qL2 <- quantile(miR_155_targets_UTR$Transcript_length, 0.5)
qL3 <- quantile(miR_155_targets_UTR$Transcript_length, 0.75)
qL4 <- quantile(miR_155_targets_UTR$Transcript_length, 1)
#function for length to name transcripts in each quartile
TLength <- function(x) { 
  if(x <=qL1) y <- "Short"
  if(x >=qL1 & x <=qL2) y <- "Med.Short"
  if(x >=qL2 & x <=qL3) y <- "Med.Long"
  if(x >=qL3) y <- "Long"
  return(y)
}
#apply above function
miR_155_targets_UTR$Length_T_q <- sapply(miR_155_targets_UTR$Transcript_length, TLength)
#ECDF plot of length and fold change
ggplot(miR_155_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_UTR, aes(x=fold, group=Length_T_q, colour=Length_T_q), stat="ecdf", size=1) +
  scale_colour_manual(name= "Transcript Length", values=cbPalette, breaks=c("Long", "Med.Long", "Med.Short", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_lengthq.tiff", width = 5, height = 3)
#KS test for fold across length bins
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_T_q=="Long"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_T_q=="Med.Long"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_T_q=="Med.Short"], miR_155_targets_UTR$fold)
ks.test(miR_155_targets_UTR$fold[miR_155_targets_UTR$Length_T_q=="Short"], miR_155_targets_UTR$fold)

#read table containing tAI values and merge with miR_155_targets
human_tAI <- read.table("human_tAI_trim.txt")
miR_155_targets_tAI <- merge(miR_155_targets, human_tAI, by.x="RefSeq_Accession.x", by.y="V1")

#define the median of tAI and make a new column for tAI above and below median
t = median(miR_155_targets_tAI$V2)
miR_155_targets_tAI$tAI_median <- ifelse(miR_155_targets_tAI$V2 <t, "Low", "High")
#KS test for fold change difference between high and low tAI messages
ks.test(miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI_median=="High"], miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI_median=="Low"])
#ECDF plot of tAI and fold change
ggplot(miR_155_targets_tAI, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_tAI, aes(x=fold, group=tAI_median, colour=tAI_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "tAI", values=cbPalette, breaks=c("High", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_tAI.tiff", width = 5, height = 3)

#define tAI quantiles
q1_32_tAI <- quantile(miR_155_targets_tAI$V2, 0.25)
q2_32_tAI <- quantile(miR_155_targets_tAI$V2, 0.5)
q3_32_tAI <- quantile(miR_155_targets_tAI$V2, 0.75)
q4_32_tAI <- quantile(miR_155_targets_tAI$V2, 1)
#function to name based on tAI
tAI <- function(x) { 
  if(x <=q1_32_tAI) y <- "Low"
  if(x >=q1_32_tAI & x <=q2_32_tAI) y <- "Med.Low"
  if(x >=q2_32_tAI & x <=q3_32_tAI) y <- "Med.High"
  if(x >=q3_32_tAI) y <- "High"
  return(y)
}
#applies above function
miR_155_targets_tAI$tAI <- sapply(miR_155_targets_tAI$V2,tAI)
#ECDF plot of fold by tAI quartiles
ggplot(miR_155_targets_tAI, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_targets_tAI, aes(x=fold, group=tAI, colour=tAI), stat="ecdf", size=1) +
  scale_colour_manual(name="tAI", values=cbPalette, breaks=c("High", "Med.High", "Med.Low", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-155_32_tAIq.tiff", width = 5, height = 3)
#KS test for fold across tAI quartiles
ks.test(miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI=="High"], miR_155_targets_tAI$fold)
ks.test(miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI=="Med.High"], miR_155_targets_tAI$fold)
ks.test(miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI=="Med.Low"], miR_155_targets_tAI$fold)
ks.test(miR_155_targets_tAI$fold[miR_155_targets_tAI$tAI=="Low"], miR_155_targets_tAI$fold)

#BELOW IS THE ANALYSIS FOR RNA FC FOR miR-155 TARGETS

#read file for miR-155 transfection RNAseq data
miR_155_RNA <- read.delim("miR155_32hr_mRNA.txt", stringsAsFactors=FALSE)
#read file for mock transfection RNAseq data
mock_mRNA <- read.delim("mRNA_mock_32hr.txt", stringsAsFactors=FALSE)
#merge miR-155 and mock fp data
miR_155_mock_RNA <- merge(miR_155_RNA, mock_mRNA, by.x="Gene_Name", by.y="Gene_Name")
#make new column for log2 fold change in fp
miR_155_mock_RNA$fold <- log2(miR_155_mock_RNA$Expression_Level_.rpkM..x/miR_155_mock_RNA$Expression_Level_.rpkM..y)
#remove NAs, Inf and -Inf from fold
miR_155_mock_RNA <- na.omit(miR_155_mock_RNA)
miR_155_mock_RNA <- subset(miR_155_mock_RNA, fold >-Inf)
miR_155_mock_RNA <- subset(miR_155_mock_RNA, fold <Inf)
#read file for mock transfection footprint data
mock_fp <- read.delim("mock_32hr_fp.txt", stringsAsFactors=FALSE)
#merge miR-155/mock fp data with mock RNAseq data
miR_155_mock_RNA_TE <- merge(miR_155_mock_RNA, mock_fp, by.x="Gene_Name", by.y="Gene_Name")
#calculate TE (fpsubmock/RNAseqsubmock) and make a new column, TE
miR_155_mock_RNA_TE$TE <- miR_155_mock_RNA_TE$Expression_Level_.rpkM./miR_155_mock_RNA_TE$Expression_Level_.rpkM..y
#remove Inf from TE
miR_155_mock_RNA_TE <- subset(miR_155_mock_RNA_TE, TE <Inf)
#read file for miR-155 target identities from Targetscan
miR_155_TS <- read.delim("mir-155_TS.txt")
#merge miR-155/mock/TE data with miR-155 target identities
miR_155_RNA_targets <- merge(miR_155_mock_RNA_TE, miR_155_TS, by.x="Gene_Name", by.y="Ortholog.of.target.gene")
miR_155_RNA_targets <- subset(miR_155_RNA_targets, TE >0)

ggplot(miR_155_RNA_targets, aes(x=log2(TE), y=fold)) +
  geom_point() + 
  stat_smooth(method=lm) + 
  ylab("RNA Fold Change (log2)") +
  xlab("TE (log2)") +
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  scale_y_continuous(limits = c(-4, 6), expand = c(0, 0)) +
  geom_text(family="Arial", fontface="italic", x = -2.5, y = 4, label = "r = 0.06866885 , p = 0.2177") +
  theme_science()
ggsave("mRNA_fold_TE_lm.tiff",  width = 5, height = 3)
cor.test(log2(miR_155_RNA_targets$TE), miR_155_RNA_targets$fold, method="pearson")



#define x as the median of TE
xr = median(miR_155_RNA_targets$TE)
#make new column to describe TE as High or Low, above or below median
miR_155_RNA_targets$TE.high.low <- ifelse(miR_155_RNA_targets$TE <xr, "Low", "High")
#KS test for fold repression
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="High"], miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="Low"])
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="High"], miR_155_RNA_targets$fold)
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="Low"], miR_155_RNA_targets$fold)
#make an ECDF plot of fold for all miR-155 targets and those with High or Low TE
ggplot(miR_155_RNA_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_RNA_targets, aes(x=fold, group=TE.high.low, colour=TE.high.low), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Low", "All")) +
  xlab("RNAseq Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-155_RNA_TE.tiff", width = 5, height = 3)
median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="High"])
median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE.high.low=="Low"])
median(miR_155_RNA_targets$fold)

#define quartile values TE
q1_TE_32 <- quantile(miR_155_RNA_targets$TE, 0.25)
q2_TE_32 <- quantile(miR_155_RNA_targets$TE, 0.5)
q3_TE_32 <- quantile(miR_155_RNA_targets$TE, 0.75)
q4_TE_32 <- quantile(miR_155_RNA_targets$TE, 1)
#function to bin TE by quartile
TE <- function(x) { 
  if(x <=q1_TE_32) y <- "Low"
  if(x >=q1_TE_32 & x <=q2_TE_32) y <- "Med.Low"
  if(x >=q2_TE_32 & x <=q3_TE_32) y <- "Med.High"
  if(x >=q3_TE_32) y <- "High"
  return(y)
}
#applies above function to new column
miR_155_RNA_targets$TE_q <- sapply(miR_155_RNA_targets$TE,TE)
#make an ECDF plot of fold for all miR-155 targets and those binned by TE quartile
ggplot(miR_155_RNA_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_155_RNA_targets, aes(x=fold, group=TE_q, colour=TE_q), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Med.High", "Med.Low", "Low", "All")) +
  xlab("RNAseq Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-155_RNA_TEq.tiff", width = 5, height = 3)
#MWW test for fold repression, TE by quartiles compared to all
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="High"], miR_155_RNA_targets$fold)
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Med.High"], miR_155_RNA_targets$fold)
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Med.Low"], miR_155_RNA_targets$fold)
ks.test(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Low"], miR_155_RNA_targets$fold)

median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="High"])
median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Med.High"])
median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Med.Low"])
median(miR_155_RNA_targets$fold[miR_155_RNA_targets$TE_q=="Low"])
median(miR_155_RNA_targets$fold)

miR_155_RNA_targets$Number_of_Sites <- miR_155_RNA_targets$Conserved.sites.total + miR_155_RNA_targets$Poorly.conserved.sites.total

ggplot(subset(miR_155_RNA_targets, Number_of_Sites <3), aes(y = fold, x = as.factor(Number_of_Sites), fill = TE.high.low)) + 
  geom_boxplot() +
  scale_fill_manual(name = "TE", values = c("#E69F00", "#56B4E9")) +
  ylab("RNAseq Fold Change (log2)") +
  xlab("Number of miR-155 Sites") +
  theme_science()
#save above plot
ggsave("miR-155_RNA_box.tiff", width = 5, height = 3)

#merge RPF and RNAseq data
translation_RNA_miR_155 <- merge(miR_155_targets, miR_155_RNA_targets, by.x="Gene_Name", by.y="Gene_Name")
#make new plot for FC RPF vs FC RNA and save plot
ggplot(translation_RNA_miR_155, aes(x=fold.y, y=fold.x)) + 
  geom_point() + 
  geom_point() + 
  xlab("Fold Change RNAseq (log2)") +
  geom_hline(yintercept = 0, size=0.7) +
  geom_vline(xintercept = 0, size=0.7) +
  ylab("Fold Change RPF (log2)") +
  stat_smooth(method=lm)  +
  geom_text(family="Arial", fontface="italic", x = -1.4, y = 6, label = "r = 0.3852321") +
  geom_text(family="Arial", fontface="italic", x = -1.4, y = 5, label = "p = 3.282e-15") +
  theme_science()
ggsave("miR_155_RNA_translation.tiff", width =4.5, height = 4.5)
#find Pearson correlation for FC RPF and FC RNA
cor.test(translation_RNA_miR_155$fold.x, translation_RNA_miR_155$fold.y, method="spearman")

#BELOW IS THE ANALYSIS OF miR-1 RPF FC

#read file for miR-1 transfection footprint data
miR_1_fp <- read.delim("mir1_32hr_fp.txt", stringsAsFactors=FALSE)
#read file for mock transfection footprint data
mock_fp <- read.delim("mock_32hr_fp.txt", stringsAsFactors=FALSE)
#merge miR-1 and mock fp data
miR_1_mock <- merge(miR_1_fp, mock_fp, by.x="Gene_Name", by.y="Gene_Name")
#make new column for log2 fold change in fp
miR_1_mock$fold <- log2(miR_1_mock$Expression_Level_.rpkM..x/miR_1_mock$Expression_Level_.rpkM..y)
#remove NAs, Inf and -Inf from fold
miR_1_mock <- na.omit(miR_1_mock)
miR_1_mock <- subset(miR_1_mock, fold >-Inf)
miR_1_mock <- subset(miR_1_mock, fold <Inf)
#read file for mock transfection RNAseq data
mock_mRNA <- read.delim("mRNA_mock_32hr.txt", stringsAsFactors=FALSE)
#merge miR-1/mock fp data with mock RNAseq data
miR_1_TE <- merge(miR_1_mock, mock_mRNA, by.x="Gene_Name", by.y="Gene_Name")
#calculate TE (fpsubmock/RNAseqsubmock) and make a new column, TE
miR_1_TE$TE <- miR_1_TE$Expression_Level_.rpkM..y/miR_1_TE$Expression_Level_.rpkM.
#remove Inf from TE
miR_1_TE <- subset(miR_1_TE, TE <Inf)
#read file for miR-1 target identities from Targetscan
miR_1_TS <- read.delim("mir-1_TS.txt")
miR_1_TS <- subset(miR_1_TS, Representative.miRNA == "hsa-miR-1-3p")
#merge miR-1/mock/TE data with miR-1 target identities
miR_1_targets <- merge(miR_1_TE, miR_1_TS, by.x="Gene_Name", by.y="Ortholog.of.target.gene")

ggplot(miR_1_targets, aes(x=log2(TE), y= fold)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  ylab("RPF Fold Change (log2)") +
  xlab("TE (log2)") +
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  scale_y_continuous(limits = c(-7, 6), expand = c(0, 0)) +
  geom_text(family="Arial", fontface="italic", x = 2.5, y = 4, label = "r = -0.3191454  , p = 1.238e-13") +
  theme_science()
ggsave("fold_logTE_lm_miR1.tiff",  width = 5, height = 3)
cor.test(log2(miR_1_targets$TE), miR_1_targets$fold, method="pearson")


miR_1_targets <- subset(miR_1_targets, TE >0)
#define x as the median of TE
x = median(miR_1_targets$TE)
#make new column to describe TE as High or Low, above or below median
miR_1_targets$TE.high.low <- ifelse(miR_1_targets$TE <x, "Low", "High")
#KS test for fold repression
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High"], miR_1_targets$fold[miR_1_targets$TE.high.low=="Low"])
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High"], miR_1_targets$fold)
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="Low"], miR_1_targets$fold)
#make an ECDF plot of fold for all miR-1 targets and those with High or Low TE
ggplot(miR_1_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets, aes(x=fold, group=TE.high.low, colour=TE.high.low), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-1_32_TE.tiff", width = 5, height = 3)

median(miR_1_targets$fold[miR_1_targets$TE.high.low=="High"])
median(miR_1_targets$fold[miR_1_targets$TE.high.low=="Low"])
median(miR_1_targets$fold)

#define quartile values TE
q1_TE_32 <- quantile(miR_1_targets$TE, 0.25)
q2_TE_32 <- quantile(miR_1_targets$TE, 0.5)
q3_TE_32 <- quantile(miR_1_targets$TE, 0.75)
q4_TE_32 <- quantile(miR_1_targets$TE, 1)
#function to bin TE by quartile
TE <- function(x) { 
  if(x <=q1_TE_32) y <- "Low"
  if(x >=q1_TE_32 & x <=q2_TE_32) y <- "Med.Low"
  if(x >=q2_TE_32 & x <=q3_TE_32) y <- "Med.High"
  if(x >=q3_TE_32) y <- "High"
  return(y)
}
#applies above function to new column
miR_1_targets$TE_q <- sapply(miR_1_targets$TE,TE)
#make an ECDF plot of fold for all miR-1 targets and those binned by TE quartile
ggplot(miR_1_targets, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets, aes(x=fold, group=TE_q, colour=TE_q), stat="ecdf", size=1) +
  scale_colour_manual(name="TE", values=cbPalette, breaks=c("High", "Med.High", "Med.Low", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save above plot
ggsave("miR-1_32_TEq.tiff", width = 5, height = 3)

#KS test for fold repression, TE by quartiles compared to all
ks.test(miR_1_targets$fold[miR_1_targets$TE_q=="High"], miR_1_targets$fold)
ks.test(miR_1_targets$fold[miR_1_targets$TE_q=="Med.High"], miR_1_targets$fold)
ks.test(miR_1_targets$fold[miR_1_targets$TE_q=="Med.Low"], miR_1_targets$fold)
ks.test(miR_1_targets$fold[miR_1_targets$TE_q=="Low"], miR_1_targets$fold)
#find the median fold change for messages with each TE quartile
median(miR_1_targets$fold[miR_1_targets$TE_q=="High"])
median(miR_1_targets$fold[miR_1_targets$TE_q=="Med.High"])
median(miR_1_targets$fold[miR_1_targets$TE_q=="Med.Low"])
median(miR_1_targets$fold[miR_1_targets$TE_q=="Low"])
median(miR_1_targets$fold)

miR_1_targets$Conserved.sites.total <- as.numeric(gsub("[*]", "", miR_1_targets$Conserved.sites.total))
miR_1_targets$Poorly.conserved.sites.total <- as.numeric(gsub("[*]", "", miR_1_targets$Poorly.conserved.sites.total))

#make a new column that sums the number of conserved and poorly conserved sites
miR_1_targets$Number_of_Sites <- miR_1_targets$Conserved.sites.total + miR_1_targets$Poorly.conserved.sites.total
#make a box plot of fold change grouped by TE and binned by number of sites
ggplot(subset(miR_1_targets, Number_of_Sites <3), aes(y = fold, x = as.factor(Number_of_Sites), fill = TE.high.low)) + 
  geom_boxplot() +
  scale_fill_manual(name = "TE", values = c("#E69F00", "#56B4E9")) +
  ylab("RPF Fold Change (log2)") +
  xlab("Number of miR-1 Sites") +
  theme_science()
#save above plot
ggsave("miR-1_32_box.tiff", width = 5, height = 3)
#KS test of the above plot
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High" & miR_1_targets$Number_of_Sites == 1], miR_1_targets$fold[miR_1_targets$TE.high.low=="Low" & miR_1_targets$Number_of_Sites == 1])
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High" & miR_1_targets$Number_of_Sites == 2], miR_1_targets$fold[miR_1_targets$TE.high.low=="Low" & miR_1_targets$Number_of_Sites == 2])
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High" & miR_1_targets$Number_of_Sites == 3], miR_1_targets$fold[miR_1_targets$TE.high.low=="Low" & miR_1_targets$Number_of_Sites == 3])
ks.test(miR_1_targets$fold[miR_1_targets$TE.high.low=="High" & miR_1_targets$Number_of_Sites == 4], miR_1_targets$fold[miR_1_targets$TE.high.low=="Low" & miR_1_targets$Number_of_Sites == 4])


#merge miR-1/mock/TE/targets with human_transcripts
miR_1_targets_UTR <- merge(miR_1_targets, human_transcripts, by.x="RefSeq_Accession.x", by.y="RefSeq.mRNA")
#deine the median of 3'UTR length and make a new column for 3'UTR length above and below median
c = median(miR_1_targets_UTR$UTR_length)
miR_1_targets_UTR$UTR_median <- ifelse(miR_1_targets_UTR$UTR_length <c, "Short", "Long")
#KS test for fold change difference between lon and short 3'UTR messages
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$UTR_median=="Long"], miR_1_targets_UTR$fold[miR_1_targets_UTR$UTR_median=="Short"])
#ECDF plot of 3'UTR length and fold change
ggplot(miR_1_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_UTR, aes(x=fold, group=UTR_median, colour=UTR_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "3'UTR Length", values=cbPalette, breaks=c("Long", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_UTR.tiff", width = 5, height = 3)
#define 3'UTR length quartiles
q1_32 <- quantile(miR_1_targets_UTR$UTR_length, 0.25)
q2_32 <- quantile(miR_1_targets_UTR$UTR_length, 0.5)
q3_32 <- quantile(miR_1_targets_UTR$UTR_length, 0.75)
q4_32 <- quantile(miR_1_targets_UTR$UTR_length, 1)
#function for 3'UTR length to name transcripts in each quartile
Length <- function(x) { 
  if(x <=q1_32) y <- "Short"
  if(x >=q1_32 & x <=q2_32) y <- "Med.Short"
  if(x >=q2_32 & x <=q3_32) y <- "Med.Long"
  if(x >=q3_32) y <- "Long"
  return(y)
}
#applies the above function
miR_1_targets_UTR$Length_q <- sapply(miR_1_targets_UTR$UTR_length, Length)
#ECDF plot of 3'UTR length and fold change
ggplot(miR_1_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_UTR, aes(x=fold, group=Length_q, colour=Length_q), stat="ecdf", size=1) +
  scale_colour_manual(name= "3'UTR Length", values=cbPalette, breaks=c("Long", "Med.Long", "Med.Short", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_UTRq.tiff", width = 5, height = 3)
#ks tests for fold change between 3'UTR length quartiles
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_q=="Long"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_q=="Med.Long"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_q=="Med.Short"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_q=="Short"], miR_1_targets_UTR$fold)

#deine the median of transcript length and make a new column for transcript length above and below median
d = median(miR_1_targets_UTR$Transcript_length)
miR_1_targets_UTR$Transcript_median <- ifelse(miR_1_targets_UTR$Transcript_length <d, "Short", "Long")
#KS test for fold change difference between long and short messages
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Transcript_median=="Long"], miR_1_targets_UTR$fold[miR_1_targets_UTR$Transcript_median=="Short"])
#ECDF plot of length and fold change
ggplot(miR_1_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_UTR, aes(x=fold, group=Transcript_median, colour=Transcript_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "Transcript Length", values=cbPalette, breaks=c("Long", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_length.tiff", width = 5, height = 3)
#define length quartiles
qL1 <- quantile(miR_1_targets_UTR$Transcript_length, 0.25)
qL2 <- quantile(miR_1_targets_UTR$Transcript_length, 0.5)
qL3 <- quantile(miR_1_targets_UTR$Transcript_length, 0.75)
qL4 <- quantile(miR_1_targets_UTR$Transcript_length, 1)
#function for length to name transcripts in each quartile
TLength <- function(x) { 
  if(x <=qL1) y <- "Short"
  if(x >=qL1 & x <=qL2) y <- "Med.Short"
  if(x >=qL2 & x <=qL3) y <- "Med.Long"
  if(x >=qL3) y <- "Long"
  return(y)
}
#apply above function
miR_1_targets_UTR$Length_T_q <- sapply(miR_1_targets_UTR$Transcript_length, TLength)
#ECDF plot of length and fold change
ggplot(miR_1_targets_UTR, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_UTR, aes(x=fold, group=Length_T_q, colour=Length_T_q), stat="ecdf", size=1) +
  scale_colour_manual(name= "Transcript Length", values=cbPalette, breaks=c("Long", "Med.Long", "Med.Short", "Short", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_lengthq.tiff", width = 5, height = 3)
#KS test for fold across length bins
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_T_q=="Long"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_T_q=="Med.Long"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_T_q=="Med.Short"], miR_1_targets_UTR$fold)
ks.test(miR_1_targets_UTR$fold[miR_1_targets_UTR$Length_T_q=="Short"], miR_1_targets_UTR$fold)

#merge tAI values with miR_1_targets
miR_1_targets_tAI <- merge(miR_1_targets, human_tAI, by.x="RefSeq_Accession.x", by.y="V1")

#deine the median of tAI and make a new column for tAI above and below median
t = median(miR_1_targets_tAI$V2)
miR_1_targets_tAI$tAI_median <- ifelse(miR_1_targets_tAI$V2 <t, "Low", "High")
#KS test for fold change difference between high and low tAI messages
ks.test(miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI_median=="High"], miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI_median=="Low"])
#ECDF plot of tAI and fold change
ggplot(miR_1_targets_tAI, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_tAI, aes(x=fold, group=tAI_median, colour=tAI_median), stat="ecdf", size=1) +
  scale_colour_manual(name= "tAI", values=cbPalette, breaks=c("High", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_tAI.tiff", width = 5, height = 3)

#define tAI quantiles
q1_32_tAI <- quantile(miR_1_targets_tAI$V2, 0.25)
q2_32_tAI <- quantile(miR_1_targets_tAI$V2, 0.5)
q3_32_tAI <- quantile(miR_1_targets_tAI$V2, 0.75)
q4_32_tAI <- quantile(miR_1_targets_tAI$V2, 1)
#function to name based on tAI
tAI <- function(x) { 
  if(x <=q1_32_tAI) y <- "Low"
  if(x >=q1_32_tAI & x <=q2_32_tAI) y <- "Med.Low"
  if(x >=q2_32_tAI & x <=q3_32_tAI) y <- "Med.High"
  if(x >=q3_32_tAI) y <- "High"
  return(y)
}
#applies above function
miR_1_targets_tAI$tAI <- sapply(miR_1_targets_tAI$V2,tAI)
#ECDF plot of fold by tAI quartiles
ggplot(miR_1_targets_tAI, aes(x=fold, colour="All")) + 
  stat_ecdf(size=1) + 
  geom_step(data=miR_1_targets_tAI, aes(x=fold, group=tAI, colour=tAI), stat="ecdf", size=1) +
  scale_colour_manual(name="tAI", values=cbPalette, breaks=c("High", "Med.High", "Med.Low", "Low", "All")) +
  xlab("RPF Fold Change (log2)") +
  ylab("Cumulative Fraction") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_science()
#save the above plot
ggsave("miR-1_32_tAIq.tiff", width = 5, height = 3)
#KS test for fold across tAI quartiles
ks.test(miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI=="High"], miR_1_targets_tAI$fold)
ks.test(miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI=="Med.High"], miR_1_targets_tAI$fold)
ks.test(miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI=="Med.Low"], miR_1_targets_tAI$fold)
ks.test(miR_1_targets_tAI$fold[miR_1_targets_tAI$tAI=="Low"], miR_1_targets_tAI$fold)




#working with UTRs
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

#get UTRs sequences and lengths
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
utr5p <- biomaRt::getSequence(id=miR_155_targets_all$RefSeq_Accession.x,type="refseq_mrna",seqType="5utr",mart=ensembl)
utr3p <- biomaRt::getSequence(id=miR_155_targets_all$RefSeq_Accession.x,type="refseq_mrna",seqType="3utr",mart=ensembl)
#put lenghts
utr3p$length3p<-nchar(utr3p$`3utr`)
utr5p$length5p<-nchar(utr5p$`5utr`)
#merge with the rest
tmp<-merge(miR_155_targets_all, utr3p, by.x="RefSeq_Accession.x", by.y="refseq_mrna")
miR_155_TE_UTR<-merge(tmp, utr5p, by.x="RefSeq_Accession.x", by.y="refseq_mrna")


#dump when sequence is unavailable

miR_155_TE_UTR <- plyr::rename(miR_155_TE_UTR, c("3utr"="utr3p", "5utr"="utr5p"))
miR_155_TE_UTR <- subset(miR_155_TE_UTR, utr5p != "Sequence unavailable")
miR_155_TE_UTR <- subset(miR_155_TE_UTR, utr3p != "Sequence unavailable")

#calculate GC

miR_155_TE_UTR$GC_utr3p<-0
miR_155_TE_UTR$GC_utr5p<-0

for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$GC_utr3p <- sum(letterFrequency(DNAString(xrow$utr3p), letters=c("G","C"), as.prob=TRUE)) }
for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$GC_utr5p <- sum(letterFrequency(DNAString(xrow$utr5p), letters=c("G","C"), as.prob=TRUE)) }


#function for calculation of MFE for each sequence
rnafold <-
  function(x){
    write(x, file="/tmp/aqq")
    fold_r=system2('RNAfold', args=c('-i','/tmp/aqq'), stdout='/tmp/out')
    fold_r=read.csv(file = "/tmp/out", header=FALSE)
    nindex=grep("[\\(\\)]",s2c(as.character(fold_r[2, ])))
    nlength=length(nindex)
    index01=nindex[nlength]-1
    index02=nindex[nlength-1]+1
    fold_f=as.numeric(paste(s2c(as.character(fold_r[2, ]))[index02:index01],collapse=""))
    return(fold_f)
  }

miR_155_TE_UTR$MFE_utr3p<-0
miR_155_TE_UTR$MFE_utr5p<-0

for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$MFE_utr5p <- rnafold(xrow$utr5p) }
for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$MFE_utr3p <- rnafold(xrow$utr3p) }

#normalization of MFE from the paper by Trotta, 2014 (PLoS One)
mfedenref<-read.csv(file="mfeden_ref.txt", header=TRUE, sep="\t")
mfeden_length<-mfedenref$L..nt.
mfeden_energy<-mfedenref$MFE..kcal.mol.
plot(mfedenref)
lm(mfeden_energy ~ mfeden_length)
abline(7.7410, -0.3117)

#got regression line, now to normalize MFE

get_MFEden<-
  function(l) {
    mfeden_normalized<- -0.3117*l + 7.7410
    return(mfeden_normalized)
  }

finalMFE<-
  function(l,mfe) {
    lzero<-8
    mfeden<-100* (mfe - get_MFEden(l))/(l - lzero)
    return(mfeden)
  }

miR_155_TE_UTR$MFE_utr3p_norm<-0
miR_155_TE_UTR$MFE_utr5p_norm<-0

for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$MFE_utr3p_norm <- finalMFE(xrow$length3p, xrow$MFE_utr3p) }
for (index in 1:nrow(miR_155_TE_UTR)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$MFE_utr5p_norm <- finalMFE(xrow$length5p, xrow$MFE_utr5p) }

miR_155_TE_UTR$Number_of_Sites <- miR_155_TE_UTR$Conserved.sites.total + miR_155_TE_UTR$Poorly.conserved.sites.total

#analysis of translational efficiency

miR_155_mRNA <- read.delim("miR155_32hr_mRNA.txt", stringsAsFactors=FALSE)
miR_155_mRNA <- plyr::rename(miR_155_mRNA, c("Expression_Level_.rpkM."="mRNA_explvl_mir155", "Read_Count"="mRNA_rc_mir155","Feature_Length"="mRNA_fl_mir155", "Total_Mapped_Reads"="mRNA_tmr_mir155"))


miR_155_TE_UTR_all <- merge(miR_155_TE_UTR, miR_155_mRNA, by.x="RefSeq_Accession.x", by.y="RefSeq_Accession")

#calculate repressed TE
miR_155_TE_UTR_all$TE_mir155 <- miR_155_TE_UTR_all$Expression_Level_.rpkM..x/miR_155_TE_UTR_all$mRNA_explvl_mir155
#remove Inf from TE
miR_155_TE_UTR_all <- subset(miR_155_TE_UTR_all, TE_mir155 <Inf)

ggplot(miR_155_TE_UTR_all, aes(x=TE, y=TE_mir155, color=fold)) + geom_point() + coord_cartesian(xlim=c(0,10), ylim=c(0,10)) + geom_smooth() + scale_colour_gradientn(colours=rainbow(4)) + theme_science()

ggplot(miR_155_TE_UTR_all, aes(x=TE, y=TE_mir155, color=fold)) + geom_point() + scale_y_continuous(limits=c(0,8)) + scale_x_continuous(limits = c(0,8))  + scale_colour_gradientn(colours=rainbow(4)) + theme_science() + geom_abline(slope = 1, intercept = 0) + geom_smooth(method = lm, se=FALSE)

miR_155_TE_UTR_all$mrnadiff<-miR_155_TE_UTR_all$Expression_Level_.rpkM.-miR_155_TE_UTR_all$mRNA_explvl_mir155
ggplot(miR_155_TE_UTR_all, aes(x=Expression_Level_.rpkM., y=mRNA_explvl_mir155, color=fold)) + 
    geom_point() + scale_colour_gradientn(colours=rainbow(5)) + 
    theme_science() + scale_y_log10() + 
    scale_x_log10() + geom_abline(slope=1, size=1) + 
    geom_abline(slope=1, intercept=0.2) + geom_abline(slope=1, intercept=-0.2)


#predicting rpfs using function derived by Eureqa
predictrpfmir <- function(rpfmock,mrnamock,mrnamir) {
  predicted<- (rpfmock*mrnamir*6.49428014187871)/(rpf_mock + 5.27420349826627*mrnamock)
  return(predicted)
}



for (index in 1:nrow(miR_155_TE_UTR_all)) { xrow = miR_155_TE_UTR[index, ]; miR_155_TE_UTR[index,]$predicted <- predictrpfmir(xrow$Expression_Level_.rpkM..y, xrow$Expression_Level_.rpkM., xrow$mRNA_explvl_mir155) }

#predicting FC using mRNA only with function derived by Eureqa
predictFC_mrna <- function(mrnamock, mrnamir) {
  predictedFCmrna <- (log2(as.numeric(mrnamir)) - log2(0.185 + as.numeric(mrnamock)))
  return(predictedFCmrna)
}


#predicting FC using mRNA and TE with function derived by Eureqa
predictFC_mrnaTE <- function(mrnamock, mrnamir, TE) {
  predictedFCmrnaTE <- (log2(mrnamir + 2.37353161434869/TE) - log2(4.6022630598024 + mrnamock + -0.207626470782861/(3.11141247407809 - mrnamock)))
  return(predictedFCmrnaTE)
}

miR_155_TE_UTR_all$predictedFC_mrna<-0
miR_155_TE_UTR_all$predictedFC_mrnaTE<-0

for (index in 1:nrow(miR_155_TE_UTR_all)) { xrow = miR_155_TE_UTR_all[index, ]; miR_155_TE_UTR_all[index,]$predictedFC_mrna <- predictFC_mrna(xrow$Expression_Level_.rpkM., xrow$mRNA_explvl_mir155) }

for (index in 1:nrow(miR_155_TE_UTR_all)) { xrow = miR_155_TE_UTR_all[index, ]; miR_155_TE_UTR_all[index,]$predictedFC_mrnaTE <- predictFC_mrnaTE(xrow$Expression_Level_.rpkM., xrow$mRNA_explvl_mir155, xrow$TE) }



ggplot(miR_155_TE_UTR_all, aes(x=fold, y= predictedFC_mrna)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  ylab("Predicted Fold Change RPF using mRNA levels") +
  xlab("Fold Change RPF (log2)") +
  geom_text(family="Arial", fontface="italic", x = 3, y = -3, label = "r = 0.4368066, p < 2.2e-16") +
  theme_science()
ggsave("FC_lm_mRNAonly.tiff",  width = 5, height = 5)

ggplot(miR_155_TE_UTR_all, aes(x=fold, y= predictedFC_mrnaTE)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  ylab("Predicted Fold Change RPF using mRNA levels and TE") +
  xlab("Fold Change RPF (log2)") +
  geom_text(family="Arial", fontface="italic", x = 3, y = -3, label = "r = 0.7474624, p < 2.2e-16") +
  theme_science()
ggsave("FC_lm_mRNATE.tiff",  width = 5, height = 5)


ggplot(miR_155_TE_UTR_all, aes(x=fold, y= MFE_utr5p_norm)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  ylab("Normalized MFE of 5' UTR") +
  xlab("Fold Change RPF") +
  geom_text(family="Arial", fontface="italic", x = 3, y = 50, label = "r = 0.006603787, p = 0.3227") +
  theme_science()
ggsave("FC_lm_MFE_5UTR.tiff",  width = 5, height = 5)

ggplot(miR_155_TE_UTR_all, aes(x=fold, y= MFE_utr3p_norm)) + 
  geom_point() + 
  stat_smooth(method=lm) + 
  geom_hline(aes(yintercept=0), size=0.7) +
  geom_vline(aes(xintercept=0), size=0.7) +
  ylab("Normalized MFE of 3' UTR") +
  xlab("Fold Change RPF") +
  geom_text(family="Arial", fontface="italic", x = 3, y = -10, label = "r = 0.02544013, p = 0.256") +
  theme_science()
ggsave("FC_lm_MFE_3UTR.tiff",  width = 5, height = 5)

cor.test(miR_155_TE_UTR_all$fold, miR_155_TE_UTR_all$predictedFC_mrna)
cor.test(miR_155_TE_UTR_all$fold, miR_155_TE_UTR_all$predictedFC_mrnaTE)


foldmrna<-lm(fold ~ predictedFC_mrna, miR_155_TE_UTR_all)

summary(foldmrna)

foldmrnaTE<-lm(fold ~ predictedFC_mrnaTE, miR_155_TE_UTR_all)

summary(foldmrnaTE)

### targeted genes vs controls

x=median(miR_155_changes_control$TE_ko_2hR)
miR_155_changes_control$TE_ko_2hR.high.low<-ifelse(miR_155_changes_control$TE_ko_2hR<x, "Low", "High")
x=median(miR_155_changes_control$TE_ko_4hR)
miR_155_changes_control$TE_ko_4hR.high.low<-ifelse(miR_155_changes_control$TE_ko_4hR<x, "Low", "High")
x=median(miR_155_changes_control$TE_ko_8hR)
miR_155_changes_control$TE_ko_8hR.high.low<-ifelse(miR_155_changes_control$TE_ko_8hR<x, "Low", "High")
x=median(miR_155_changes_control$TE_ko_48hR)
miR_155_changes_control$TE_ko_48hR.high.low<-ifelse(miR_155_changes_control$TE_ko_48hR<x, "Low", "High")

ggplot(miR_155_changes_targets, aes(x=TE_ko_2hR*RPF.2.h, group=TE_2hr_.high.low, color=TE_2hr_.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Targets at 2hrs")
ggsave("S14_targets_2hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_targets, aes(x=TE_ko_4hR*RPF.4.h, group=TE_4hr_.high.low, color=TE_4hr_.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4))  + labs(title="Targets at 4hrs")
ggsave("S14_targets_4hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_targets, aes(x=TE_ko_8hR*RPF.8.h, group=TE_8hr_.high.low, color=TE_8hr_.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Targets at 8hrs")
ggsave("S14_targets_8hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_targets, aes(x=TE_ko_48hR*RPF.48.h, group=TE_48hr_.high.low, color=TE_48hr_.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Targets at 48hrs")
ggsave("S14_targets_48hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_control, aes(x=TE_ko_2hR*RPF.2.h, group=TE_ko_2hR.high.low, color=TE_ko_2hR.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Control at 2hrs")
ggsave("S14_control_2hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_control, aes(x=TE_ko_4hR*RPF.4.h, group=TE_ko_4hR.high.low, color=TE_ko_4hR.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Control at 4hrs")
ggsave("S14_control_4hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_control, aes(x=TE_ko_8hR*RPF.8.h, group=TE_ko_8hR.high.low, color=TE_ko_8hR.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Control at 8hrs")
ggsave("S14_control_8hrs.tiff",  width = 5, height = 5)

ggplot(miR_155_changes_control, aes(x=TE_ko_48hR*RPF.48.h, group=TE_ko_48hR.high.low, color=TE_ko_48hR.high.low)) + 
  geom_density(alpha=.3) + theme_science() + 
  scale_x_continuous(limits = c(-4,4)) + labs(title="Control at 48hrs")
ggsave("S14_control_48hrs.tiff",  width = 5, height = 5)

target<-miR_155_changes_targets$TE_ko_2hR*miR_155_changes_targets$RPF.2.h
control<-miR_155_changes_control$TE_ko_2hR*miR_155_changes_control$RPF.2.h
ks.test(target, control, alternative = "greater")
target<-miR_155_changes_targets$TE_ko_4hR*miR_155_changes_targets$RPF.4.h
control<-miR_155_changes_control$TE_ko_4hR*miR_155_changes_control$RPF.4.h
ks.test(target, control, alternative = "greater")
target<-miR_155_changes_targets$TE_ko_8hR*miR_155_changes_targets$RPF.8.h
control<-miR_155_changes_control$TE_ko_8hR*miR_155_changes_control$RPF.8.h
ks.test(target, control, alternative = "greater")
target<-miR_155_changes_targets$TE_ko_48hR*miR_155_changes_targets$RPF.48.h
control<-miR_155_changes_control$TE_ko_48hR*miR_155_changes_control$RPF.48.h
ks.test(target, control, alternative = "greater")

### tAI analysis - testing some weakly documented features



#trna raw files contain the numbers for respective tRNAs and zeros for stop codons (last three lines)


#this is the order from codon file; 
#column  codon
#-------------
#1	TTT
#2	TTC
#3	TTA
#4	TTG
#...
#first case, TTT codon corresponds to TTT tRNA (wrong)
dros.trna<-scan("../code/codon_analysis/drosophila.trna.raw.trna")


#this is the correct order, where tRNAs are complemented: TTT codon corresponds to AAA tRNA, TTC codon corresponds to GAA tRNA, etc.
dros.compl<-scan("../code/codon_analysis/drosophila.trna.raw.complement")

#let's calculate both ways

dros.ws.trna<-get.ws(tRNA = dros.trna, sking=0) #sking set to zero, indicating eukaryota
dros.ws.compl<-get.ws(tRNA = dros.compl, sking=0) #sking set to zero, indicating eukaryota

#all CDS from Ensembl
dros.m<-matrix(scan("../code/codon_analysis/drosophila2.m"), ncol=61, byrow = TRUE) 
dros.m<-dros.m[,-33] 
dros.tai <- get.tai(dros.m, dros.ws.trna)
dros.tai.compl <- get.tai(dros.m, dros.ws.compl)

hist(dros.tai, breaks=50)
hist(dros.tai.compl, breaks=50)
theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

dme_tAI <- as.data.frame(dros.tai.compl)

bp <- ggplot(dme_tAI, aes(x=dros.tai.compl, fill= ..x..)) + 
  geom_histogram(aes(y=..count../sum(..count..))) + 
  scale_fill_gradient(low = "light grey", high = "black") +
  scale_y_continuous(limits = c(0, .25), expand = c(0, 0)) +
  xlab("tAI") +
  ylab("Frequency") +
  geom_vline(xintercept=c(0.298, 0.387, 0.494, 0.602), size=1) +
  theme_science()
bp + coord_flip() + 
  scale_x_reverse() + 
  theme(axis.title.y = element_text(angle = 0)) + 
  guides(fill=FALSE)
ggsave("dme_tAI_hist.tiff", height=3.1, width=3.6, units = "in")

human.trna <- scan("../code/codon_analysis/human.trna_raw.complement")
human.ws<-get.ws(tRNA = human.trna, sking=0) #sking set to zero, indicating eukaryota
human.m <- matrix(scan("../code/codon_analysis/human.m"), ncol=61, byrow = TRUE)
human.m <- human.m[,-33]
human.tai <- get.tai(human.m, human.ws)
summary(human.tai)
id.list <- scan("../code/codon_analysis/id.list",what = "character")
human.tai.values<-data.frame(id.list, human.tai)
colnames(human.tai.values)<- c("V1", "V2")
write.csv(human.tai.values, file="human_tai_values_new.csv")

ren.m<-matrix(scan("/Users/Kyle/Desktop/codonR/ren.m"), ncol=61, byrow = TRUE) 
ren.m<-ren.m[,-33] 
ren.tai.compl <- get.tai(ren.m, dros.ws.compl)


