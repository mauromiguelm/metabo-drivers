## This code analyse metabolomics data

# load packages and definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\plate_converter.R")
library(openxlsx)
library(dplyr)
library(rhdf5)

#import drug sensitity metrics

setwd(path_data_file)

data_GR50 <- read.csv("outcomes_growth_inhibition50.csv")
GR24_outliers_high <- read.csv('GR24_outliers_high.csv')
GR24_outliers_low <- read.csv('GR24_outliers_low.csv')
  
#import metabolomics 

setwd(path_metabolomics_in)

dataContent<- h5ls("metabolomics_raw.h5")

data<- rhdf5::h5read(file = "metabolomics_raw.h5", '/data')

ions <- rhdf5::h5read(file = "metabolomics_raw.h5", '/annotation')

ions <- data.frame(ions)

setwd(paste(path_data_file,'metabolomics', sep = "//"))

#import cleaned metadata

metadata <- read.csv("metadata_clean.csv")

# impot R/S groups

setwd(path_data_file)

RS_groups <- read.csv("outcomes_GR24_RSgroups_filtered.csv"  )

#import pathway defition

setwd("C:\\Users\\masierom\\polybox\\Programing\\GSEA\\GSEA_2005_updatedMauro_20180914")

source('GSEA.1.0.Mauro_20180914.R')

load('DataSource.RData')  # This should be used as mock data.


#define drugs from controls

drugs_in_screen <- c(unique(metadata$drug)[!(unique(metadata$drug) %in% c("PBS", "DMSO"))])

# calculating log2(FCs) for all data and for every drug.  --------

#iterating over plates, cells, drug, and conc to calculate FCs

metadata$cell_plate <- paste(metadata$cell, metadata$source_plate, sep = "_")

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))


lapply(unique(metadata$cell_plate), function(cell_plate_idx){
  #subset metadata to  drug & control groups
  #cell_plate_idx <- unique(metadata$cell_plate)[1]
  if (!any(grepl(cell_plate_idx,x=list.files()))){
    print(cell_plate_idx)
    results_list <- list()
    
    tmp <- subset(metadata, cell_plate == cell_plate_idx, select = c("idx","drug", "cell", "conc", "source_plate",'quadrant')) #keep one time point with GR24.. all time points have the same value..
    
    tmp_control <- subset(tmp, drug == "DMSO")
    
    print(paste("dimention of control vector =",dim(tmp_control)))
    
    tmp_drug <- subset(tmp, drug %in% drugs_in_screen)
    
    tmp_drug$drug_conc <- paste(tmp_drug$drug, tmp_drug$conc, sep = "_")
    
    #for each drug and concentration, calculate fold changes to control
    
    drug_conc <- unique(paste(tmp_drug$drug, tmp_drug$conc, sep = "_"))
    #drug_conc_idx <- "Irinotecan 0.11"
    lapply(drug_conc, function(drug_conc_idx){
      #iterate over each drug_conc
      #drug_conc_idx <- drug_conc[1]
      tmp_drug_conc <- subset(tmp_drug,drug_conc == drug_conc_idx)
      
      if(nrow(tmp_drug_conc)>0){
        
        lapply(1:nrow(data), function(metab_idx){
          #metab_idx = 1
          
          drug_metab <- data.frame("intensity"= data[metab_idx,tmp_drug_conc$idx])
          
          drug_metab$Group.1 <- tmp_drug_conc$quadrant
          
          control_metab <- aggregate(data[metab_idx,tmp_control$idx], by = list(tmp_control$quadrant),median)
          
          drug_fcs <- merge(drug_metab, control_metab, by = 'Group.1')
          
          drug_fcs <- log2(drug_fcs[,2]/ drug_fcs[,3])
          
          p_value <- t.test(drug_fcs, mu =0)$p.value
          
          return(list(paste(cell_plate_idx,drug_conc_idx, metab_idx, sep ="_"), median(drug_fcs), p_value))
        })-> results
        
        return(results)
        
      }
      
      
    }) -> results
    
    metab_fcs <- unlist(results, recursive = FALSE)
    metab_fcs <- do.call(rbind, metab_fcs)
    
    metab_fcs <- data.frame(metab_fcs)
    
    metab_fcs<- cbind(data.frame(stringr::str_split_fixed(string = metab_fcs[,1],"_",n = 5)),metab_fcs[,2:3])
    
    metab_fcs[,6] <- as.numeric(metab_fcs[,6])
    metab_fcs[,7] <- as.numeric(metab_fcs[,7])
    
    write.csv(metab_fcs, paste0(cell_plate_idx,'.csv'))
    
  }else{
    NA
  }
})


# # use linear regression to see if any metabolites have associati --------
# lm(metab~GR50) across all cell lines and concentrations that we have filtered strong effects/unnefective concentrations
# do this wiht methotrexate as we have clear expectations

# import log2fc data 

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

lapply(unique(paste(metab_fcs$drug,metab_fcs$ionIndex, sep = "_")),function(drug_ion){
  #for every ion in each drug associate with GR50
  print(drug_ion)
  #drug_ion <-  "UK5099_1"
  
  tmp <- subset(metab_fcs, drug == strsplit(drug_ion, split = "_")[[1]][1] &
                ionIndex == as.numeric(strsplit(drug_ion, split = "_")[[1]][2]))
  
  tmp_growth_metrics <- subset(data_GR50, agent == strsplit(drug_ion, split = "_")[[1]][1], select=c("cell_line", 'GR50'))
  
  if(nrow(tmp_growth_metrics) > 0 & sum(is.finite(tmp_growth_metrics$GR50))>=3){
    #only include drugs where n of defined G50 is higher or equal than 3
    
    #remove Inf and replace by a very high concentration, defining these cells as resistant
    tmp_growth_metrics$GR50  <- ifelse(tmp_growth_metrics$GR50 == Inf, max(tmp_growth_metrics$GR50[is.finite(tmp_growth_metrics$GR50)],na.rm = T)*10, tmp_growth_metrics$GR50)
    
    #combine metabolomics with GR50
    
    tmp <- merge(tmp, tmp_growth_metrics, by = "cell_line")
    
    tmp <- tmp%>% dplyr::group_by(cell_line,drug)%>% dplyr::arrange(concentration) %>% dplyr::mutate(concentration = sequence(n()))
    
    #filter too strong concentrations and ineffective concentrations
    #remove drug_conc with no effect
    
    drug_conc_to_remove <- subset(GR24_outliers_low,outliers == "low" & Drug == strsplit(drug_ion, split = "_")[[1]][1], select =c('Drug', 'Final_conc_uM') )
    
    tmp <- tmp[!paste(tmp$drug,tmp$concentration, sep = "_") %in% paste(drug_conc_to_remove$Drug,drug_conc_to_remove$Final_conc_uM, sep = "_"),]
    
    # remove cell_drug_conc with too strong effect at GR24
    
    drug_conc_to_remove <- subset(GR24_outliers_high,outliers == 'high' & Drug == strsplit(drug_ion, split = "_")[[1]][1], select =c('Drug','cell', 'Final_conc_uM') )
    
    tmp <- tmp <- tmp[!paste(tmp$cell_line,tmp$drug,tmp$concentration, sep = "_") %in% paste(drug_conc_to_remove$cell,drug_conc_to_remove$Drug,drug_conc_to_remove$Final_conc_uM, sep = "_"),]
    
    #calculate linear regression intensity vs. GI50
    
    slope <- lm(log2fc~log10(GR50), tmp)
    
    pvalue <- summary(slope)
    
    return(c(drug_ion,slope$coefficients[2],pvalue$r.squared,pvalue$adj.r.squared,pvalue$coefficients[2,4]))
    
    
  }else{
    return(NA)
  }
  
}) -> slope_metabolite_effect_on_growth


tmp <- do.call(rbind, slope_metabolite_effect_on_growth)

tmp <-data.frame(tmp)


tmp$log10.GR50. <- as.numeric(tmp$log10.GR50.)
tmp$V4 <- as.numeric(tmp$V4)
tmp$V5 <- as.numeric(tmp$V5)

#plot most interesting cases for each drug

#select the top 5 abs(log2fc) for each drug

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}


tmp3 <- stringr::str_split_fixed(string = tmp[,1],"_",n = 2)

tmp3 <- cbind(tmp3, tmp[,-c(1)])

tmp3 <- data.frame(tmp3)

colnames(tmp3) <- c('drug', 'ionIndex', 'slope', 'r2', 'adj-r2', 'pvalue')

tmp3 <- subset(tmp3, abs(slope) >= 0.08)

tmp3 <- tmp3 %>% group_by(drug) %>% dplyr::arrange(abs(slope)) %>% dplyr::slice_tail(n=5)

tmp3 <- tmp3 %>% group_by(drug) %>%dplyr::arrange(abs(slope)) %>%  dplyr::mutate(slope_norm = c(5:1)[1:n()])

tmp3$color <- ifelse(tmp3$slope<0, "Negative", "Positive")

tmp_ions <- ions[,c("ionIndex",'mzLabel','score', "name")] %>% dplyr::group_by(ionIndex) %>% dplyr::arrange(score) %>% dplyr::slice(n())

tmp3 <- merge(tmp3, tmp_ions, by='ionIndex')

tmp3$drugion <- paste(tmp3$drug, tmp3$name)

tmp3 <- tmp3[order(tmp3$drug),]


tmp3$angle_1 <- seq(0, 360, length.out = nrow(tmp3))

tmp3$radius_1 <- 1

library(ggrepel)

ggplot(tmp3, aes(x=angle_1, y=radius_1, size = abs(slope)))+
  geom_jitter(width = 2, height = 0)+
  scale_x_continuous(limits = c(0, 360))+
  scale_y_continuous(limits =c(-5,1))+
  coord_polar()+
  theme_bw()+
  geom_text(aes(label=name),size =2)+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(axis.text.x = element_text(
    angle= -90 - 360 / length(unique(tmp3$drugion)) * seq_along(tmp3$drugion))) -> p


p+ geom_bar(aes(x=angle_1, y=-1, col = factor(drug)),width = 2, height = 0, size = 3, stat = 'identity')+
  geom_text(aes(x=angle_1, y=0,label=drug),size =2)+
  coord_polar()

#plot with circlize

df <- tmp3
library(circlize)

df$x=0
df$color <- ifelse(df$slope<0, "#FF0000", "#00FF00")
circos.par("track.height" = 0.3,cell.padding = c(0.02, 0.04, 0.02, 0.04))
setwd(path_fig)

par(mar = c(500, 800, 400, 200) + 100)
png("association_GI50_metabolite.png",width = 10000,height = 10000,res = 700) 

circos.initialize(df$drug, x = (df$slope_norm))

circos.track(df$drug, y =df$slope,
             panel.fun = function(x, y) {
               circos.axis(major.tick = F,labels = NULL)
  },bg.border = NA)

circos.trackText(sectors = df$drug,x = df$slope_norm, y= df$x, labels = df$name, 
                 cex= 0.8,track.index = 1,facing = 'clockwise',col = df$color,adj =  c(0),niceFacing = T)

circos.track(df$drug, y =df$slope,
             panel.fun = function(x, y) {
               circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)


circos.trackPoints(sectors = df$drug,x = df$slope_norm, y= df$x+2,cex = abs(df$slope)+0.5, col = df$color,pch = 16) 


circos.track(df$drug, y =df$slope,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(1), 
                           CELL_META$sector.index,cex = 0.01,facing = 'clockwise',niceFacing = T)
             },bg.border = NA)


color = viridis::viridis(nrow(df))
i=1
for(x in (df$drug)){
  highlight.sector(track.index = 2,col = color[i],sector.index = x)
  highlight.sector(track.index = 3, col= color[i], sector.index = x)
  i=i+1
  
  
}

circos.track(df$drug, y =df$slope,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(1)+3, 
                           CELL_META$sector.index,cex = 0.9,facing = 'clockwise',niceFacing = T)
             },bg.border = NA)

dev.off()
circos.clear()
#plot iwth ggplot2

df <- tmp3[,c('drug','name','slope','angle_1')]
library(tidyverse)

df <- df %>%dplyr::group_by(drug) %>% dplyr::mutate(value = n())

colnames(df) <- c("name", "type",'slope','angle_1','value')

df <- df %>% dplyr::filter(type != "all") %>%
  dplyr::mutate(name = as.factor(name))%>%
  dplyr::arrange(name, type) %>%
  dplyr::mutate(type = as.factor(type)) #fct_reorder2(name, value))

lvl0 <- tibble(name = "Parent", value = 0, fill = NA, fill2=NA, angle_1 = 0,angle_2 = NA,level = 0)

lvl1 <- df %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(value = sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fill = name)%>%
  dplyr::mutate(fill2 = NA)%>%
  dplyr::mutate(angle_1 = seq(90, 450, length.out = n()))%>%
  dplyr::mutate(level = 1)
  
lvl2 <- df %>%
  dplyr::select(name = type, value, fill = name, fill2 =slope) %>%
  dplyr::mutate(level = 2)

lvl2$angle_1 <- seq(90, 450, length.out = nrow(lvl2))

library(ggnewscale)

df <- bind_rows(lvl0, lvl1, lvl2)

df$name <- as.character(df$name)
df$fill <- as.character(df$fill)
df$fill <- as.factor(df$fill)

df %>%
  #mutate(name = as.factor(name) %>% fct_reorder2(fill, value)) %>%
  #arrange(as.factor(name)) %>%
  mutate(level = as.factor(level)) %>%
  ggplot(aes(x = level, y = value, fill = fill, alpha = level)) +
  geom_col(width = 1, color = "gray90", size = 0.25, position = position_stack()) +
  #geom_text(aes(label = name, angle =), size = 2, position = position_stack(vjust = 0.5)) +
  scale_colour_brewer(fill, type = 'div')+
  new_scale("fill")+
  geom_col(aes(x = level, y = value, fill = fill2),width = 1, color = "gray90", size = 0.25, position = position_stack(),inherit.aes = T) +
  scale_fill_gradient2("fill2",low = "blue", mid = "white", high = "red")+
  coord_polar(theta = "y",start = 0) +
  geom_text(aes(label = name, angle = angle_1), size = 2, position = position_stack(vjust = 0.5),hjust=0.5) +
  #scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 0.7), guide = F) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  #scale_fill_brewer(palette = "Dark2", na.translate = F) +
  labs(x = NULL, y = NULL) +
  theme_minimal()


library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

# Create dataset
data <- data.frame(
  individual=tmp3$drug,
  value=tmp3$name
)

# Set a number of 'empty bar'
empty_bar <- 10

# Add lines to the initial dataset
to_add <- matrix(NA, empty_bar, ncol(data))
colnames(to_add) <- colnames(data)
data <- rbind(data, to_add)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make the plot
p <- ggplot(data, aes(x=as.factor(value),y=id, label=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_text() +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) + 
  geom_text(data=label_data, aes(x=id, y=value, label=individual, hjust=hjust), color="black", fontface="bold", size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p



# plot one interesting case

doi <- "Panzem-2-ME2"
iot <- 103

tmp2 <- subset(metab_fcs, ionIndex == iot & drug == doi)

tmp2$experiment <- paste(tmp2$cell_line, tmp2$drug)

tmp2 <- merge(tmp2, data_GR50[,c('experiment', "GR50")], by = "experiment")

tmp2$GR50  <- ifelse(tmp2$GR50 == Inf, max(tmp2$GR50[is.finite(tmp2$GR50)],na.rm = T)*10, tmp2$GR50)

#filter too strong concentrations and ineffective concentrations
#remove drug_conc with no effect

tmp2 <- tmp2%>% dplyr::group_by(cell_line,drug)%>% dplyr::arrange(concentration) %>% dplyr::mutate(concentration = sequence(n()))

drug_conc_to_remove <- subset(GR24_outliers_low,outliers == "low" & Drug == doi, select =c('Drug', 'Final_conc_uM') )

tmp2 <- tmp2[!paste(tmp2$drug,tmp2$concentration, sep = "_") %in% paste(drug_conc_to_remove$Drug,drug_conc_to_remove$Final_conc_uM, sep = "_"),]

# remove cell_drug_conc with too strong effect at GR24

drug_conc_to_remove <- subset(GR24_outliers_high,outliers == 'high' & Drug == doi, select =c('Drug','cell', 'Final_conc_uM') )

tmp2 <- tmp2 <- tmp2[!paste(tmp2$cell_line,tmp2$drug,tmp2$concentration, sep = "_") %in% paste(drug_conc_to_remove$cell,drug_conc_to_remove$Drug,drug_conc_to_remove$Final_conc_uM, sep = "_"),]

ggplot(tmp2, aes(y=log2fc,x=log10(GR50), label = cell_line))+
  geom_point(col = 'red')+
  geom_text_repel()

# RS groups ---------------------------------------------------------------


# #calculate pathway enrichment for R/S groups ------------------------------------

lapply(unique(RS_groups$Drug), function(drug_idx){
  #drug_idx = 'BPTES'
  print(drug_idx)
  RS_sub <- subset(RS_groups, Drug == drug_idx & !is.na(percent_change_GR))
  
  res_sub <- subset(RS_sub, group == 'R' )
    
  sens_sub <-subset(RS_sub, group == 'S')
  
  meta_sub <- subset(metadata, drug == drug_idx)
  
  #TODO apply this code to the previous definition of conc_idxs
  meta_sub$conc <- meta_sub %>% dplyr::group_by(conc) %>%  group_indices(conc)
  #unique(meta_sub$tmp)
  #xtabs(~conc+tmp, meta_sub)
  
  #TODO add "I" to balance cases
  resistant_metadata_idx <- subset(meta_sub, cell %in% unique(res_sub$cell) & conc %in%  unique(res_sub$Final_conc_uM))
  
  sensitive_metadata_idx <- subset(meta_sub, cell %in% unique(sens_sub$cell) & conc %in%  unique(sens_sub$Final_conc_uM))
  
  max_length <- min(nrow(resistant_metadata_idx),nrow(sensitive_metadata_idx))
  
  if(max_length>6){
    sensitive_metadata_idx <- sensitive_metadata_idx[1:max_length,]
    
    resistant_metadata_idx <- resistant_metadata_idx[1:max_length,]
    
    class.v <- rep(0,nrow(sensitive_metadata_idx))
    
    dataset <- data[,sensitive_metadata_idx$idx]
    
    dataset <- cbind(dataset, data[,resistant_metadata_idx$idx])
    
    dataset <- data.frame(dataset)
    
    ions_sub <- ions[,c('ionIndex','idKEGG')]
    
    ions_sub <- tidyr::separate_rows(ions_sub,idKEGG, convert = TRUE, sep = ' ')
    
    ions_sub <- subset(ions_sub,!grepl("^\\s*$", idKEGG))
    
    ions_sub <- ions_sub %>% group_by(idKEGG) %>% slice(1)
    
    dataset <- dataset[ions_sub$ionIndex,]
    
    rownames(dataset) <- ions_sub$idKEGG
    class.v <- append(class.v,rep(1,nrow(resistant_metadata_idx)))
    
    CLS <- list(class.v = class.v, phen = c("S", 'R'))
    
    
    setwd("C:\\Users\\masierom\\polybox\\Programing\\GSEA\\GSEA_2005_updatedMauro_20180914")
    
    Output_GSEA <- GSEA(
      dataset = dataset,                       # Input gene expression Affy dataset file in RES or GCT format
      CLS = CLS,                               # Input class vector (phenotype) file in CLS format
      temp =  temp,                            # Gene set database in GMT format
      output.directory      = getwd(),         # Directory where to store output and results (default: "")
      #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
      doc.string            = drug_idx,      # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
      non.interactive.run   = T,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
      reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
      nperm                 = 500,            # Number of random permutations (default: 1000)
      weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
      nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
      fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
      fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
      topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
      adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
      gs.size.threshold.min = 10,               # Minimum size (in genes) for database gene sets to be considered (default: 25)
      gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
      reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
      preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
      random.seed           = 760435,          # Random number generator seed. (default: 123456)
      perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
      fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
      replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
      save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
      OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
      use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
    )
      
  }
  
  return(Output_GSEA)
  
}) -> Output_GSEA

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')



# Compare metabolomics results for R/Sgroups, and see which one is better

# plot the metabolites that are increased in R or S as a volcano plot



# check which metabolites are driving the correlation, hopefully positive controls
# Correlate FC with GR50 and GR24
# For metabolic drugs: check if direct substrate/products involved in MoA are relating to any drug metric 
# Multiomics: Link protein/gene/mRNA levels to metabolomics/drug sensitivity
# Multiomics: distance_to_target analysis