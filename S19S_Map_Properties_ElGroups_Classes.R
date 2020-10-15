#### Plotting ICR sample data on a map
# Kings Creek (N04D) in Kansas had a latitude listed as ~69...which is not in Kansas. The other King's Creek lat ~39, so we
# adjusted accordingly

# Load in libraries
library(reshape2) # For reordering data
library(ggplot2) # For pretty plots
library(ggpubr) # Combining plots into a single figure
library(stringr) # String manipulation
library(dplyr) # Miscellaneous modifications

# ################## #
#### Load in data ####
# ################## #

setwd("path/to/Metabolome Data Files")

# Load in ICR data
data = read.csv("Processed_S19S_Sed-Water_08.12_Data.csv", row.names = 1, stringsAsFactors = F)
mol = read.csv("Processed_S19S_Sed-Water_08.12_Mol.csv", row.names = 1, stringsAsFactors = F)

# Load in meta-data
S19S.lat.long = read.csv("S19S_Metadata.csv", stringsAsFactors = F)

# Ordering meta data
S19S.lat.long = S19S.lat.long[order(S19S.lat.long$ID),]

# Subsetting only to lat.long
S19S.lat.long = S19S.lat.long[,c("ID", "Latitude", "Longitude")]

# Renaming first column
colnames(S19S.lat.long)[1] = "Site"

# Selecting only molecular formulae
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Calculating StoC ratios
mol$StoC_ratio = with(mol, S/C)

# Fixing compound class definitions
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"


# ################################################################ #
#### Combining samples from the same site for plotting purposes ####
# ################################################################ #

# Finding surface and sediment samples
sed.samples = colnames(data)[grep("Field", colnames(data))]
wat.samples = colnames(data)[-grep("Field|INC", colnames(data))]

# Finding unique samples
sed.samples = unique(gsub("_ICR.*$", "", sed.samples))
wat.samples = unique(gsub("_ICR.[0-9].*$", "_ICR", wat.samples))

unique.samples = c(wat.samples, sed.samples)

rm(sed.samples, wat.samples)

# Creating new matrix
derep.data = matrix(nrow = nrow(data), ncol = length(unique.samples), dimnames = list(row.names(data), unique.samples))
derep.data = as.data.frame(derep.data)

for(i in 1:ncol(derep.data)){
  temp = as.matrix(data[,grep(colnames(derep.data)[i], colnames(data))])
  
  if(ncol(temp) > 1){
    derep.data[,i] = rowSums(temp)
  } else {
    derep.data[,i] = temp
  }
}

# Creating factors sheet
derep.factors = data.frame(Samples = colnames(derep.data), Type = "Surface Water", stringsAsFactors = F) # Adding sample type
derep.factors$Type[grep("Field", derep.factors$Samples)] = "Sediment"
derep.factors$Site = str_extract(derep.factors$Samples, "S19S_[0-9]{4}") # Adding sample site
derep.factors = derep.factors %>% left_join(S19S.lat.long, by = "Site") # Adding in lat/longs


# ####################################################### #
#### Calculating average derived statistics by sample #####
# ####################################################### #

char = data.frame(NOSC = NA, DBE = NA, AI_Mod = NA, NtoC = NA, PtoC = NA, StoC = NA,
                  Peaks = NA, Type = derep.factors$Type, Site = derep.factors$Site, 
                  Latitude = derep.factors$Latitude, Longitude = derep.factors$Longitude,
                  row.names = colnames(derep.data), stringsAsFactors = F)

for(i in 1:ncol(derep.data)){
  temp = mol[which(derep.data[,i] > 0),]
  
  char$NOSC[i] = mean(temp$NOSC, na.rm = T)
  char$DBE[i] = mean(temp$DBE, na.rm = T)
  char$AI_Mod[i] = mean(temp$AI_Mod, na.rm = T)
  char$NtoC[i] = mean(temp$NtoC_ratio, na.rm = T)
  char$PtoC[i] = mean(temp$PtoC_ratio, na.rm = T)
  char$StoC[i] = mean(temp$StoC_ratio, na.rm = T)
  char$Peaks[i] = nrow(temp)
  
  rm(temp)
} # Looping through samples and averaging characteristics

rm(i)

# Separating surface and sediment samples
surf = char[which(char$Type %in% "Surface Water"),]
sed = char[which(char$Type %in% "Sediment"),]

# Contiguous US boundaries
north = 49.3457868
south = 24.7433195
east = -66.9513812
west = -124.7844079

missi = -90.19789 # Longitude where the Mississippi River cross through St. Louis, MI, USA

# Filtering data based on US boundaries
sed.us = sed[which(sed$Longitude < east & sed$Longitude > west & sed$Latitude > south & sed$Latitude < north),]
surf.us = surf[which(surf$Longitude < east & surf$Longitude > west & surf$Latitude > south & surf$Latitude < north),]


### Performing East/West analyses
# Idenifying east-west of the Mississippi River
sed.us$East.West = "East"
sed.us$East.West[which(sed.us$Longitude < missi)] = "West"

surf.us$East.West = "East"
surf.us$East.West[which(surf.us$Longitude < missi)] = "West"

# Melting and combining data
melt.var = rbind(melt(sed.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")),
                 melt(surf.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")))

# Plotting east-west comparisons
melt.var = melt.var[which(melt.var$variable %in% c("NOSC", "DBE", "AI_Mod", "Peaks")),]
ggplot(melt.var, aes(x = East.West, y = value))+
  geom_boxplot()+facet_grid(variable~Type, scales = "free_y")+
  xlab(NULL)+theme_bw()+theme(axis.text = element_text(color = "black", size = 12),
                              axis.text.x = element_text(color = "black", size = 14),
                              axis.ticks = element_line(color = "black"),
                              panel.border = element_rect(color = "black"),
                              panel.background = element_blank(),
                              panel.grid = element_blank())

# East/West characteristic stats
stats = NULL

for(curr.var in unique(as.character(melt.var$variable))){
  # Sediment stats
  w = which(colnames(sed.us) %in% curr.var)
  temp = data.frame(Variable = sed.us[,w], East.West = sed.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Sediment East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water stats
  w = which(colnames(surf.us) %in% curr.var)
  temp = data.frame(Variable = surf.us[,w], East.West = surf.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Water East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water vs. Sediment stats
  w = which(colnames(char) %in% curr.var)
  temp = data.frame(Variable = char[,w], Type = char$Type)
  temp = wilcox.test(Variable~Type, data = temp)
  temp = data.frame(Comparison = paste("Water vs. Sed.", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
}

rm(curr.var, w, temp)

### Generating maps
for(var in c("NOSC", "DBE", "AI_Mod", "NtoC", "PtoC", "StoC", "Peaks")){
  # Selecting desired characteristic
  sed = data.frame(Longitude = sed.us$Longitude, Latitude = sed.us$Latitude, 
                   Variable = sed.us[,which(colnames(sed.us) %in% var)])
  surf = data.frame(Longitude = surf.us$Longitude, Latitude = surf.us$Latitude, 
                    Variable = surf.us[,which(colnames(surf.us) %in% var)])
  
  # Limits on ggplot scale
  limits = c(min(c(sed$Variable, surf$Variable)), 
             max(c(sed$Variable, surf$Variable)))
  
  # Plotting US data
  sed.plot = ggplot(data = sed)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = Variable), size = 7, alpha = 0.7)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "dodgerblue2", mid = "goldenrod2",
                          high = "firebrick2", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Sediment", var)) + labs(color = var)+
    labs(color = paste0("Average ", var))+
    coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.text.y=element_blank(),
                                       axis.ticks=element_blank(),
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       panel.background=element_blank(),
                                       panel.border=element_blank(),
                                       panel.grid.major=element_blank(),
                                       panel.grid.minor=element_blank(),
                                       plot.background=element_blank())
  
  surf.plot = ggplot(data = surf)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = Variable), size = 7, alpha = 0.7)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "dodgerblue2", mid = "goldenrod2",
                          high = "firebrick2", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Surface Water", var)) + labs(color = var)+
    labs(color = paste0("Average ", var))+
    coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.text.y=element_blank(),
                                       axis.ticks=element_blank(),
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       panel.background=element_blank(),
                                       panel.border=element_blank(),
                                       panel.grid.major=element_blank(),
                                       panel.grid.minor=element_blank(),
                                       plot.background=element_blank())
  
  ggarrange(sed.plot, surf.plot, ncol = 1, common.legend = T) %>%
    ggexport(filename = paste0("Maps 08-21-20/US_Maps_Metric_", var, ".pdf"))
}


# ################################################# #
#### Calculating elemental composition by sample ####
# ################################################# #

el.comp = matrix(data = 0, nrow = ncol(derep.data), ncol = length(unique(mol$El_comp)), 
                 dimnames = list(colnames(derep.data), unique(mol$El_comp)))

for(i in 1:nrow(el.comp)){
  temp = mol[which(derep.data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(el.comp)){
    el.comp[i,j] = length(which(temp$El_comp %in% colnames(el.comp)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

el.comp = as.data.frame(t(apply(el.comp, 1, function(x) (x/sum(x))*100)))
el.comp$Site = derep.factors$Site
el.comp$Longitude = derep.factors$Longitude
el.comp$Latitude = derep.factors$Latitude
el.comp$Type = "Surface Water"
el.comp$Type[grep("Field", row.names(el.comp))] = "Sediment"

surf.comp = el.comp[-grep("Field", row.names(el.comp)),]
sed.comp = el.comp[grep("Field", row.names(el.comp)),]

# Filtering data based on US boundaries
sed_comp.us = sed.comp[which(sed.comp$Longitude < east & sed.comp$Longitude > west & sed.comp$Latitude > south & sed.comp$Latitude < north),]
surf_comp.us = surf.comp[which(surf.comp$Longitude < east & surf.comp$Longitude > west & surf.comp$Latitude > south & surf.comp$Latitude < north),]

### Performing East/West analyses
# Idenifying east-west of the Mississippi River
sed_comp.us$East.West = "East"
sed_comp.us$East.West[which(sed_comp.us$Longitude < missi)] = "West"

surf_comp.us$East.West = "East"
surf_comp.us$East.West[which(surf_comp.us$Longitude < missi)] = "West"

# Melting and combining data
melt.var = rbind(melt(sed_comp.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")),
                 melt(surf_comp.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")))
melt.var$variable = factor(melt.var$variable, levels = c("CHO", "CHON", "CHOS", "CHOP", "CHONS", "CHONP",
                                                         "CHOSP", "CHONSP"))

# Plotting east-west comparisons
ggplot(melt.var, aes(x = East.West, y = value))+
  geom_boxplot()+facet_grid(variable~Type, scales = "free_y")+
  xlab(NULL)+theme_bw()+theme(axis.text = element_text(color = "black", size = 12),
                              axis.text.x = element_text(color = "black", size = 14),
                              axis.ticks = element_line(color = "black"),
                              panel.border = element_rect(color = "black"),
                              panel.background = element_blank(),
                              panel.grid = element_blank())

# Calculating East/West stats
for(curr.var in unique(as.character(melt.var$variable))){
  # Sediment east/west stats
  w = which(colnames(sed_comp.us) %in% curr.var)
  temp = data.frame(Variable = sed_comp.us[,w], East.West = sed_comp.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Sediment East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water east/west stats
  w = which(colnames(surf_comp.us) %in% curr.var)
  temp = data.frame(Variable = surf_comp.us[,w], East.West = surf_comp.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Water East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water vs. Sediment stats
  w = which(colnames(el.comp) %in% curr.var)
  temp = data.frame(Variable = el.comp[,w], Type = el.comp$Type)
  temp = wilcox.test(Variable~Type, data = temp)
  temp = data.frame(Comparison = paste("Water vs. Sed.", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
}

rm(curr.var, w, temp)

### Generating Maps
for(el.var in c("CHO", "CHON", "CHONS", "CHONP", "CHONSP", "CHOP", "CHOS", "CHOSP")){
  # Selecting desired elemental composition
  sed = data.frame(Longitude = sed_comp.us$Longitude, Latitude = sed_comp.us$Latitude, 
                   ElComp = sed_comp.us[,which(colnames(sed_comp.us) %in% el.var)])
  surf = data.frame(Longitude = surf_comp.us$Longitude, Latitude = surf_comp.us$Latitude, 
                    ElComp = surf_comp.us[,which(colnames(surf_comp.us) %in% el.var)])
  
  # Limits on ggplot scale
  limits = c(min(c(sed$ElComp, surf$ElComp)), 
             max(c(sed$ElComp, surf$ElComp)))
  
  # Plot US data
  sed.plot = ggplot(data = sed)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = ElComp), size = 7, alpha = 0.7)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "darkgreen", mid = "goldenrod2",
                          high = "purple", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Sediment", el.var, "Occurence")) + labs(color = el.var)+
    labs(color = paste0("%", el.var))+
    coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.text.y=element_blank(),
                                       axis.ticks=element_blank(),
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       panel.background=element_blank(),
                                       panel.border=element_blank(),
                                       panel.grid.major=element_blank(),
                                       panel.grid.minor=element_blank(),
                                       plot.background=element_blank())
  
  surf.plot = ggplot(data = surf)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = ElComp), size = 7, alpha = 0.7)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "darkgreen", mid = "goldenrod2",
                          high = "purple", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Surface Water", el.var, "Occurence")) + labs(color = el.var)+
    labs(color = paste0("%", el.var))+
    coord_fixed() +  theme_bw() + theme(axis.line=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.text.y=element_blank(),
                                        axis.ticks=element_blank(),
                                        axis.title.x=element_blank(),
                                        axis.title.y=element_blank(),
                                        panel.background=element_blank(),
                                        panel.border=element_blank(),
                                        panel.grid.major=element_blank(),
                                        panel.grid.minor=element_blank(),
                                        plot.background=element_blank())
  
  ggarrange(sed.plot, surf.plot, ncol = 1, common.legend = T) %>%
    ggexport(filename = paste0("Maps 08-21-20/US_Maps_El_", el.var, ".pdf"))
}


# #################################################### #
#### Calculating compound classifications by sample ####
# #################################################### #

comp.class = matrix(data = 0, nrow = ncol(derep.data), ncol = length(unique(mol$bs1_class)), 
                 dimnames = list(colnames(derep.data), unique(mol$bs1_class)))

for(i in 1:nrow(comp.class)){
  temp = mol[which(derep.data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(comp.class)){
    comp.class[i,j] = length(which(temp$bs1_class %in% colnames(comp.class)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

comp.class = as.data.frame(t(apply(comp.class, 1, function(x) (x/sum(x))*100)))
comp.class$Site = derep.factors$Site
comp.class$Longitude = derep.factors$Longitude
comp.class$Latitude = derep.factors$Latitude
comp.class$Type = "Surface Water"
comp.class$Type[grep("Field", row.names(el.comp))] = "Sediment"

surf.class = comp.class[-grep("Field", row.names(comp.class)),]
sed.class = comp.class[grep("Field", row.names(comp.class)),]

# Filtering data based on US boundaries
sed_class.us = sed.class[which(sed.class$Longitude < east & sed.class$Longitude > west & sed.class$Latitude > south & sed.class$Latitude < north),]
surf_class.us = surf.class[which(surf.class$Longitude < east & surf.class$Longitude > west & surf.class$Latitude > south & surf.class$Latitude < north),]

### Performing East/West analyses
# Idenifying east-west of the Mississippi River
sed_class.us$East.West = "East"
sed_class.us$East.West[which(sed_class.us$Longitude < missi)] = "West"

surf_class.us$East.West = "East"
surf_class.us$East.West[which(surf_class.us$Longitude < missi)] = "West"

# Melting and combining data
melt.var = rbind(melt(sed_class.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")),
                 melt(surf_class.us, id.vars = c("Type", "Site", "Latitude", "Longitude", "East.West")))

# Plotting east-west comparisons
ggplot(melt.var, aes(x = East.West, y = value))+
  geom_boxplot()+facet_grid(variable~Type, scales = "free_y")+
  xlab(NULL)+theme_bw()+theme(axis.text = element_text(color = "black", size = 12),
                              axis.text.x = element_text(color = "black", size = 14),
                              axis.ticks = element_line(color = "black"),
                              panel.border = element_rect(color = "black"),
                              panel.background = element_blank(),
                              panel.grid = element_blank())

# Calculating East/West stats
for(curr.var in unique(as.character(melt.var$variable))){
  # Sediment east/west stats
  w = which(colnames(sed_class.us) %in% curr.var)
  temp = data.frame(Variable = sed_class.us[,w], East.West = sed_class.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Sediment East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water east/west stats
  w = which(colnames(surf_class.us) %in% curr.var)
  temp = data.frame(Variable = surf_class.us[,w], East.West = surf_class.us$East.West)
  temp = wilcox.test(Variable~East.West, data = temp)
  temp = data.frame(Comparison = paste("Water East/West", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
  
  # Water vs. Sediment stats
  w = which(colnames(comp.class) %in% curr.var)
  temp = data.frame(Variable = comp.class[,w], Type = comp.class$Type)
  temp = wilcox.test(Variable~Type, data = temp)
  temp = data.frame(Comparison = paste("Water vs. Sed.", curr.var), W = temp$statistic, p.value = temp$p.value)
  stats = rbind(stats, temp)
}

stats$p.value = p.adjust(stats$p.value, method = "fdr")
stats$p.value = round(stats$p.value, digits = 5)
stats = stats[-which(stats$p.value > 0.05),]

rm(curr.var, w, temp)

### Generating Maps
for(class.var in c("Amino Sugar", "Carbohydrate", "Cond Hydrocarbon", "Lignin", "Lipid", "Protein", "Tannin", "Unsat Hydrocarbon")){
  # Selecting desired elemental composition
  sed = data.frame(Longitude = sed_class.us$Longitude, Latitude = sed_class.us$Latitude, 
                   Class = sed_class.us[,which(colnames(sed_class.us) %in% class.var)])
  surf = data.frame(Longitude = surf_class.us$Longitude, Latitude = surf_class.us$Latitude, 
                    Class = surf_class.us[,which(colnames(surf_class.us) %in% class.var)])
  
  # Limits on ggplot scale
  limits = c(min(c(sed$Class, surf$Class)), 
             max(c(sed$Class, surf$Class)))
  
  # Plot US data
  sed.plot = ggplot(data = sed)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = Class), size = 7, alpha = 0.8)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "goldenrod", mid = "gray",
                          high = "cyan4", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Sediment", class.var, "Occurence")) + labs(color = class.var)+
    labs(color = paste0("%", class.var))+
    coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.text.y=element_blank(),
                                       axis.ticks=element_blank(),
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       panel.background=element_blank(),
                                       panel.border=element_blank(),
                                       panel.grid.major=element_blank(),
                                       panel.grid.minor=element_blank(),
                                       plot.background=element_blank())
  
  surf.plot = ggplot(data = surf)+ 
    borders("state", colour = "black", fill = "white")+ 
    geom_point(aes(x = Longitude, y = Latitude, color = Class), size = 7, alpha = 0.8)+
    geom_vline(xintercept = missi, color = "red", lty = 2, lwd = 2)+
    scale_color_gradient2(limits = limits, low = "goldenrod", mid = "gray",
                          high = "cyan4", midpoint = (max(limits)+min(limits))/2)+
    ggtitle(paste("Surface Water", class.var, "Occurence")) + labs(color = class.var)+
    labs(color = paste0("%", class.var))+
    coord_fixed() +  theme_bw() + theme(axis.line=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.text.y=element_blank(),
                                        axis.ticks=element_blank(),
                                        axis.title.x=element_blank(),
                                        axis.title.y=element_blank(),
                                        panel.background=element_blank(),
                                        panel.border=element_blank(),
                                        panel.grid.major=element_blank(),
                                        panel.grid.minor=element_blank(),
                                        plot.background=element_blank())
  
  ggarrange(sed.plot, surf.plot, ncol = 1, common.legend = T) %>%
    ggexport(filename = paste0("Maps/US_Maps_Class_", class.var, ".pdf"))
}

rm(sed, surf)