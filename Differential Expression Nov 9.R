#Tri's code from before (thanks!)
library(dplyr)
library(ggplot2)
library(data.table)


df = readRDS("C:\\Users\\Faith's PC\\Documents\\Data 467\\Project\\GSE39612_edat.RDS")
df2 <- data.frame(t(df))
row_names = row.names(df2)
df2$Patient = row_names
rownames(df2) = NULL

y = read.csv("C:\\Users\\Faith's PC\\Documents\\Data 467\\Project\\patient_data.csv", header = FALSE, col.names = c('Patient','Tumor_type', 'Sample_type', 'Sample_name'))
y$Tumor_type = gsub(',','', y$Tumor_type)
y$Sample_type = gsub(',','', y$Sample_type)
y$Sample_name = gsub(',','', y$Sample_name)

merged_df = merge(df2,y, by = 'Patient')
numerical_df = merged_df[,1:20184]
categorical_df = merged_df[,c(1,20185,20186,20187)]
merged_df; numerical_df; categorical_df


#Necessary packages
library(limma)
library(edgeR)

#Getting just the group of samples that were renormalized together (no replicates)
renorm_df = df[, -(1:40)]
renorm_cat_df = categorical_df[-(1:40),]


#Differential expression for MCC (cancer) and normal skin (non-cancerous)
#This section is for the linear model for the data frame
group = factor(renorm_cat_df$Tumor_type)
design = model.matrix(~0 + group)
colnames(design) = levels(group)
fit = lmFit(renorm_df, design)

#This section gets the comparison that we're actually interested in, 
#then uses empirical Bayes to reduce variance based on the rest of the data frame
contr = makeContrasts(MCC - normal, levels = colnames(fit))
tmp = contrasts.fit(fit, contr)
tmp = eBayes(tmp)

#This section helps us arrange the table into something more readable
table = topTable(tmp, sort.by = "P", n = Inf)
table$Gene = rownames(table)
table = table[,c("Gene", names(table)[1:6])]
write.table(table, file = "C:\\Users\\Faith's PC\\Documents\\Data 467\\Project\\DE MCC vs Normal.txt", sep = "\t")


#Differential expression for primary (less severe) and metastic (more severe).
group = factor(renorm_cat_df$Sample_type)
design = model.matrix(~0 + group)
colnames(design) = levels(group)
fit = lmFit(renorm_df, design)

contr = makeContrasts(metastatic - primary, levels = colnames(fit))
tmp = contrasts.fit(fit, contr)
tmp = eBayes(tmp)

table = topTable(tmp, sort.by = "P", coef = 1, n = Inf)
table$Gene = rownames(table)
table = table[,c("Gene", names(table)[1:6])]
write.table(table, file = "C:\\Users\\Faith's PC\\Documents\\Data 467\\Project\\DE Primary vs Metastatic.txt", sep = "\t")

