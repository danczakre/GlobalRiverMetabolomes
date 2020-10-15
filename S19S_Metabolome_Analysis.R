### Demonstration R script to plot whole S19S FTICR-metric distributions

# Swtiches
match.wat.sed = F # Controls whether or not the sites are matched

# Load in necessary libraries first
library(ggplot2); library(reshape2) # For pretty plots
library(vegan) # For ecology
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation

# ######################### #
#### Data Pre-processing ####
# ######################### #

# Set working directory
setwd("/path/to/Metabolome Data Files")

# Load in mol. data (contains values we want to plot)
data = read.csv("Processed_S19S_Sed-Water_08.12_Data.csv", stringsAsFactors = F, row.names = 1)
mol = read.csv("Processed_S19S_Sed-Water_08.12_Mol.csv", stringsAsFactors = F, row.names = 1)

# Removing poorly calibrated samples from the dataset
poor.cal = read.csv(list.files(pattern = "*_Poorly_Calibrated_Samples.csv"))

if(length(poor.cal[,1]) > 1){
  data = data[,-which(colnames(data) %in% gsub("-", ".", poor.cal$samples))]
} else {
  stop("You've specified to remove poor calibrants, but there were none provided.")
}

rm(poor.cal)

# Error checking step that ensures the masses in my mol file match my data file exactly
if(!identical(row.names(mol), row.names(data))){
  stop("Your masses in the data file and mol file do not match. Maybe incorrect files were loaded.")
}

# Removing peaks that were not assigned a molecular formula
na.loc = which(mol$MolForm %in% NA)
mol = mol[-na.loc,]
data = data[-na.loc,]
rm(na.loc)

# Given that ICR data cannot reliable track abundance/concentration, we need to set data to presence/absence
data[data > 0] = 1

# Renaming bs1_classes
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"

# Adding mass
mol$Mass = as.numeric(as.character(row.names(mol)))

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Sample_Type = "Surface Water", stringsAsFactors = F)
factors$Sample_Type[grep("Field", factors$Samples)] = "Sediment"
factors$Site = str_extract(factors$Samples, "[0-9]{4}")

### Matching sites, if the switch is set
if(match.wat.sed == T){
  
  # Matching data across samples types
  surf.sites = unique(factors$Site[which(factors$Sample_Type %in% "Surface Water")]) # Sites in surface water
  sed.sites = unique(factors$Site[which(factors$Sample_Type %in% "Sediment")]) # Sites in sediment
  common.sites = intersect(surf.sites, sed.sites) # Sites common to both sample types
  
  # Subset factors based upon common sites
  factors = factors[which(factors$Site %in% common.sites),]
  data = data[,which(colnames(data) %in% factors$Samples)]
  
  # Drop missing molecular formula
  mol = mol[-which(rowSums(data) == 0),]
  data = data[-which(rowSums(data) == 0),]
}

rm(match.wat.sed)

# ##################################################################### #
#### Analyzing all formula found in either surface water or sediment ####
# ##################################################################### #

# Unique sample types (i.e., water and sediment)
uniq.sample.type = unique(factors$Sample_Type)

# Shorthand way of doing those 6 commands above
data.by.type = as.data.frame(matrix(nrow = nrow(data), ncol = length(uniq.sample.type), 
                                    dimnames = list(row.names(data), uniq.sample.type)))

# Using a for-loop to run through our data and identify peaks present in each sample type
for(i in 1:length(uniq.sample.type)){
  # Find temporary datasets for sample types
  temp = data[,which(factors$Sample_Type %in% uniq.sample.type[i])]
  
  # Add in sums into the parent data frame
  data.by.type[,i] = rowSums(temp)
}

rm(temp, i)

# Resetting data back to presence/absence
data.by.type[data.by.type > 0] = 1

# Partitioning the molecular information based upon sample
temp.mol = mol[which(data.by.type$`Surface Water` > 0),] # Finding molecular information for surface water
temp.mol = temp.mol[,c("AI_Mod", "DBE", "NOSC")] # Selecting subset of interesting variables
temp.mol = melt(as.matrix(temp.mol)) # Melting as a matrix to get the Var1/Var2 melt format
temp.mol$Sample_Type = "Surface Water" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = temp.mol # Creating the object to eventually go into ggplot

temp.mol = mol[which(data.by.type$`Sediment` > 0),]
temp.mol = temp.mol[,c("AI_Mod", "DBE", "NOSC")] # Selecting subset of interesting variables
temp.mol = melt(as.matrix(temp.mol)) # Melting as a matrix to get the Var1/Var2 melt format
temp.mol$Sample_Type = "Sediment" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = rbind(melt.mol.by.type, temp.mol)

rm(temp.mol)

# Adding in molecular formula-by-sample count into our melt.mol.by.type object
melt.mol.by.type = rbind(melt.mol.by.type, data.frame(Var1 = colnames(data), 
                                                      Var2 = "Molecular Formula Count",
                                                      value = colSums(data), 
                                                      Sample_Type = factors$Sample_Type))

# Statistics
metric.stats = NULL
uniq.met = unique(as.character(melt.mol.by.type$Var2))

for(i in 1:length(uniq.met)){
  w = which(melt.mol.by.type$Var2 %in% uniq.met[i])
  wil.test = wilcox.test(value~Sample_Type, data = melt.mol.by.type[w,], alternative = "two.sided")
  wil.test = data.frame(Comparison = uniq.met[i], W = wil.test$statistic, P.value = wil.test$p.value,
                        Test = wil.test$alternative)
  metric.stats = rbind(metric.stats, wil.test)
}

metric.stats$P.value = p.adjust(metric.stats$P.value, method = "fdr")

rm(wil.test, w, i)

# Plotting metrics
ggplot(melt.mol.by.type, aes(x = value, group = Sample_Type))+
  geom_density(aes(fill = Sample_Type), alpha = 0.5)+
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=12, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())


# ################################################ #
#### Comparing el. comp. and bs1 between groups ####
# ################################################ #

### Elemental Composition by sample
el.comp = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$El_comp)), 
                  dimnames = list(colnames(data), unique(mol$El_comp)))

for(i in 1:nrow(el.comp)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(el.comp)){
    el.comp[i,j] = length(which(temp$El_comp %in% colnames(el.comp)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

el.comp = as.data.frame(t(apply(el.comp, 1, function(x) (x/sum(x))*100))) # Relative abundance
el.comp = cbind(factors, el.comp)
el.comp = melt(el.comp, id.vars = colnames(factors))

### Compound Class by sample
comp.class = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$bs1_class)), 
                  dimnames = list(colnames(data), unique(mol$bs1_class)))

for(i in 1:nrow(comp.class)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(comp.class)){
    comp.class[i,j] = length(which(temp$bs1_class %in% colnames(comp.class)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

comp.class = as.data.frame(t(apply(comp.class, 1, function(x) (x/sum(x))*100)))
comp.class = cbind(factors, comp.class)
comp.class = melt(comp.class, id.vars = colnames(factors))
comp.class$variable = gsub("Hydrocarbon", "HC", comp.class$variable) # Shortening names
comp.class$variable = gsub("Carbohydrate", "Carb.", comp.class$variable)

# Performing stats on el. comp
el.stats = NULL

for(curr.el in unique(el.comp$variable)){
  temp = el.comp[which(el.comp$variable %in% curr.el),]
  max.val = max(temp$value)
  temp = wilcox.test(value~Sample_Type, data = temp, alternative = "two.sided")
  temp = data.frame(Comparison = curr.el, W = temp$statistic, p.value = temp$p.value, max.val = max.val)
  el.stats = rbind(el.stats, temp)
}

el.stats$p.value = p.adjust(el.stats$p.value, method = "fdr")
el.stats$Symbol = symnum(el.stats$p.value, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

rm(curr.el, temp)

# Performing stats on comp. class
comp.stats = NULL

for(curr.comp in unique(comp.class$variable)){
  temp = comp.class[which(comp.class$variable %in% curr.comp),]
  max.val = max(temp$value)
  temp = wilcox.test(value~Sample_Type, data = temp, alternative = "two.sided")
  temp = data.frame(Comparison = curr.comp, W = temp$statistic, p.value = temp$p.value, max.val = max.val)
  comp.stats = rbind(comp.stats, temp)
}

comp.stats$p.value = p.adjust(comp.stats$p.value, method = "fdr")
comp.stats$Symbol = symnum(comp.stats$p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

rm(curr.comp, temp)

# Making boxplots for compound class and elem. comp.
el.plot = ggplot(data = el.comp)+
  geom_boxplot(aes(x = variable, y = value, color = Sample_Type))+ 
  theme_bw() + xlab(NULL) + ylab("Relative Abundance (%)") + labs(color = "Sample Type:")+
  geom_text(data = el.stats, aes(x = Comparison, y = max.val+2, label = as.character(Symbol)), size = 7)+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

comp.plot = ggplot(data = comp.class)+
  geom_boxplot(aes(x = variable, y = value, color = Sample_Type))+ 
  theme_bw() + xlab(NULL) + ylab("Relative Abundance (%)") + labs(color = "Sample Type:")+
  geom_text(data = comp.stats, aes(x = Comparison, y = max.val+2, label = as.character(Symbol)), size = 7)+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

ggarrange(el.plot, comp.plot, ncol = 1, common.legend = T)

# ####################################################### #
#### Generating pie charts for elem. composition & bs1 ####
# ####################################################### #

# Partitioning the molecular information based upon sample
surfacew_mol = mol[which(data.by.type$`Surface Water` > 0),] # Finding molecular information for surface water
sediment_mol = mol[which(data.by.type$`Sediment` > 0),]

# Generate tables for the frequency of elemental compositions
surfacew_table_El_comp = as.data.frame(table(surfacew_mol$El_comp))
sediment_table_El_comp = as.data.frame(table(sediment_mol$El_comp))

surfacew_table_El_comp$Freq = (surfacew_table_El_comp$Freq/sum(surfacew_table_El_comp$Freq))*100
sediment_table_El_comp$Freq = (sediment_table_El_comp$Freq/sum(sediment_table_El_comp$Freq))*100

# Generate tables for the freqency of boundary sets
surfacew_table_bs1_class = as.data.frame(table(surfacew_mol$bs1_class))
sediment_table_bs1_class = as.data.frame(table(sediment_mol$bs1_class))

surfacew_table_bs1_class$Freq = (surfacew_table_bs1_class$Freq/sum(surfacew_table_bs1_class$Freq))*100
sediment_table_bs1_class$Freq = (sediment_table_bs1_class$Freq/sum(sediment_table_bs1_class$Freq))*100

# Adding sample type qualifiers 
surfacew_table_El_comp$SampleType = "Surface Water"
sediment_table_El_comp$SampleType = "Sediment"

surfacew_table_bs1_class$SampleType = "Surface Water"
sediment_table_bs1_class$SampleType = "Sediment"

# Combining frequency tables
CombinedTable_El_comp = rbind(surfacew_table_El_comp, sediment_table_El_comp)
CombinedTable_bs1_class = rbind(surfacew_table_bs1_class, sediment_table_bs1_class)

rm(list = ls(pattern = c("sediment|surfacew")))

# Calculating position
CombinedTable_El_comp = CombinedTable_El_comp %>% group_by(SampleType) %>%
  mutate(Position = (cumsum(Freq) - Freq/2)) %>% mutate(Label = paste0(round(Freq, digits = 2), "%"))
CombinedTable_bs1_class = CombinedTable_bs1_class %>% group_by(SampleType) %>%
  mutate(Position = (cumsum(Freq) - Freq/2)) %>% mutate(Label = paste0(round(Freq, digits = 2), "%"))

# Create piechart from table data
col = colorRampPalette(c("#b2182b", "white", "#2166ac"))(8)
Comp_pie = ggplot(data = CombinedTable_El_comp, aes(x = "", y = Freq, fill = Var1))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(SampleType~., ncol = 1)+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = col)+theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks =  element_blank(), panel.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank())

col = colorRampPalette(c("#762a83", "white", "#1b7837"))(10)
Class_pie = ggplot(data = CombinedTable_bs1_class, aes(x = "", y = Freq, fill = Var1))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  facet_wrap(SampleType~., ncol = 1)+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = col)+theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks =  element_blank(), panel.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank())

ggarrange(Comp_pie, Class_pie)


# ########################### #
#### Multivariate Analyses ####
# ########################### #

# Principal component analysis
pca = prcomp(x = t(data))

# Everything below this line should not change very much whether we use PCA or NMDS
ordination.scores = scores(pca) # Works with both PCA and NMDS, change the object accordingly
ordination.scores = as.data.frame(ordination.scores) # ggplot doesn't like matrices - needs to be converted to a data frame
ordination.scores$SampleType = factors$Sample_Type # Adding in sample type to our ordination scores object

# We have everything necessary for ggplot - we want to plot PC1 and PC2
ggplot(data = ordination.scores, aes(x = PC1, y = PC2, color = SampleType))+
  xlab(paste0("PC1 (", summary(pca)$importance[2,1]*100, "%)"))+
  ylab(paste0("PC2 (", summary(pca)$importance[2,2]*100, "%)"))+
  geom_point(size = 2) + theme_bw()+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

# PERMANOVA
dist = vegdist(t(data), method = "euclidean", binary = T)
perm = adonis(dist~factors$Sample_Type, permutations = 999)

# Beta-dispersion analysis
beta.disp = betadisper(dist, group = factors$Sample_Type)
beta.disp = data.frame(Type = as.character(beta.disp$group), Distance = as.numeric(beta.disp$distances),
                       stringsAsFactors = F)
beta.stats = wilcox.test(Distance~Type, data = beta.disp)

# Plotting beta-diserpsion
ggplot(data = beta.disp, aes(x = Type, y = Distance))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+ylab("Distance to Centroid")+
  theme_bw()+theme(legend.position = "none",
                   axis.text = element_text(color = "black", size = 11),
                   axis.title = element_text(color = "black", size = 13),
                   axis.ticks = element_line(color = "black"),
                   panel.border = element_rect(color = "black"),
                   panel.background = element_blank(),
                   panel.grid = element_blank())