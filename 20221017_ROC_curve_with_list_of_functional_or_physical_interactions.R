# This script is to draw ROC curve for determining cutoff of number of cell line 

setwd("C:/Rdata")
dir.create("20221017_ROC_curve_for_cutoff_with_functional_interactions")
setwd("C:/Rdata/20221017_ROC_curve_for_cutoff_with_functional_interactions")

# activate package "Epi" for drawing ROC curve
library(Epi)

# inport list of known functional interactions
# this list is located at "Dropbox/Okamura Lab share folder/Hirota/results_and_matterials/20220616_ROC_curve_for_determing_cutoff_of_cell_number"
func.list <-read.table("list_of_RBP_CCLE_miRNA_functional_interaction.txt",sep="\t",header = T,stringsAsFactors = F)
func.list[,1] <-paste0(func.list[,1],"_",func.list[,3],"_",func.list[,4],"_vs_",func.list[,5])

# inport list of Treiber's physical interactions
# this list is located at "https://github.com/Ryosuke-Hirota/20221017_ROC_curve_with_list_of_functional_or_physical_interactions"
phy.list <-read.table("list_of_treiber_physical_interaction_between_RBP_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)
phy.list[,1] <-paste0(phy.list[,1],"_",phy.list[,2],"_",phy.list[,3],"_vs_",phy.list[,4])
phy.list <-phy.list[phy.list[,5]>3,]

# inport result of correlation analysis
# this result is located at "https://github.com/Ryosuke-Hirota/20221017_ROC_curve_with_list_of_functional_or_physical_interactions"
setwd("C:/Rdata/ROC_curve_for_cutoff")
result <-read.table("combine_results_of_correlation_between_residual_and_RBP_exp.txt",sep="\t",header = T,stringsAsFactors = F)
result[,1 ] <-paste0(result[,1],"_",result[,2],"_",result[,3],"_vs_",result[,4])

# find combinations which are known to functional interactions from correlation analysis result
m <-match(func.list[,1],result[,1])
m <-na.omit(m)
result[m,10] <-1
result[-m,10] <-0

# draw ROC curve based on known functional interaction
pdf("ROC_curve_functional_interaction.pdf")

ROC(test=result$number_of_cell, stat=result$significant, plot="ROC")

dev.off()

# find combinations with physical interactions from correlation analysis result
m1 <-match(phy.list[,1],result[,1])
m1 <-na.omit(m1)
result[m1,10] <-1
result[-m1,10] <-0

# draw ROC curve based on known functional interaction
pdf("ROC_curve_physical_interaction.pdf")

ROC(test=result$number_of_cell, stat=result$significant, plot="ROC")

dev.off()

