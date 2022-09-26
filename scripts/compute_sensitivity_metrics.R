## This code takes growth curve from large screen and calculate drug metrics, 

# load packages and definitions -------------------------------------------

library(tidyr); library(ggplot2);library(gplots);library(RColorBrewer);library(viridis); library(GRmetrics)

path_metabolomics_in <- "./metetabolomics"
path_data_file = "./data"
path_fig = "./figures"
path_metadata = "./metadata"

setwd(path_data_file)
load("growth_curves_filtered.RData")

time_treatment <- 0 #when was the drug added
cutoff_GR_max_growth_effect <- 50 #pairs of drug & conc with abs change < this number will be excluded
cutoff_GR_conc_min_CV <- 10 #pairs of drug & conc with CV < this number will be excluded
drugs_in_screen <- c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO"))]))

#import log2fc metabolomics data

setwd(paste0(".",path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

#import cell metadata from depmap 

setwd(paste0("../../.",path_metadata))

cell_metadata <- read.csv("metadata_cells_depmap.csv")

#import cleaned metadata

setwd(paste0(".",path_metadata,'//','metabolomics'))
metadata <- read.csv("metadata.csv")

# calculate growth inhibition 50 metrics -----------------------------------

#split analysis between plate 1 and plate 2

lapply(c("P1",'P2'),function(plate){
  
  data_GRmetrics <- subset(data_corrected, source_plate == plate)
  
  data_Grmetrics_Ttm <- data_GRmetrics
  data_Grmetrics_Ctr <- data_GRmetrics
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO" & Final_conc_uM %in% c(367)) # I separated the drug from the control, as its easier to parse
  
  # keep only matching time points
  
  match_time_intervals <- intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time) #make sure both datasets cover the same time
  
  time_treatment <- max(match_time_intervals[match_time_intervals < time_treatment])
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)
  
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time == 72 | Time == 0)
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time == 72 | Time == 0)
  
  data_Grmetrics_Ctr$Time <- data_Grmetrics_Ctr %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
  data_Grmetrics_Ttm$Time <- data_Grmetrics_Ttm %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
  
  # Prepare data_Grmetrics_Ctr as example Case C
  
  data_Grmetrics_Ctr$time = data_Grmetrics_Ctr$Time
  
  tmp <- strsplit(x = rownames(data_Grmetrics_Ctr), "_")
  
  tmp <- lapply(tmp, "[[",1)
  
  tmp <- gsub(tmp, pattern = "/", replacement = "_")
  
  tmp <- strsplit(x = tmp, "_")
  
  tmp <- lapply(tmp, "[[",3)
  
  data_Grmetrics_Ctr$cell_line <- tmp
  
  data_Grmetrics_Ctr$agent = "-"
  
  data_Grmetrics_Ctr$concentration = 0
  
  data_Grmetrics_Ctr$cell_count <- data_Grmetrics_Ctr$Conf
  
  data_Grmetrics_Ctr <- data_Grmetrics_Ctr[,c("cell_line", "agent","time", "concentration", "cell_count")]
  
  rownames(data_Grmetrics_Ctr) <- NULL
  
  # Prepare data_Grmetrics_Ttm as example Case C
  
  head(data_Grmetrics_Ttm)
  
  data_Grmetrics_Ttm$time = data_Grmetrics_Ttm$Time
  
  tmp <- strsplit(x = rownames(data_Grmetrics_Ttm), "_")
  
  tmp <- lapply(tmp, "[[",1)
  
  tmp <- gsub(tmp, pattern = "/", replacement = "_")
  
  tmp <- strsplit(x = tmp, "_")
  
  tmp <- lapply(tmp, "[[",3)
  
  data_Grmetrics_Ttm$cell_line <- tmp; rm(tmp)
  
  data_Grmetrics_Ttm$agent <- (data_Grmetrics_Ttm$Drug)
  
  data_Grmetrics_Ttm$concentration = data_Grmetrics_Ttm$Final_conc_uM
  
  data_Grmetrics_Ttm$cell_count <- data_Grmetrics_Ttm$Conf
  
  data_Grmetrics_Ttm <- data_Grmetrics_Ttm[,c("cell_line", "agent","time", "concentration", "cell_count")]
  
  rownames(data_Grmetrics_Ttm) <- NULL
  
  # merge drug df and control df
  
  data_comb <- rbind(data_Grmetrics_Ctr, data_Grmetrics_Ttm)
  
  data_comb$cell_line <- as.character(data_comb$cell_line)

  output1 = GRmetrics::GRfit(inputData = data_comb, groupingVariables =
                               c("cell_line", 'agent'), case = "C")
  
  output<-GRmetrics::GRgetMetrics(output1)
  
  return(output)
  
  
}) -> output_IC50

output_IC50 <- do.call(rbind, output_IC50)

#save output

setwd(path_data_file)

write.csv(output_IC50, "outcomes_growth_inhibition50.csv")

# #calculate percent inhibition for GR24 ----------------------------------

lapply(c("P1", "P2"), function(plate){ 
  # calculate results by plate
  #plate = "P2"
  
  tmp <- subset(data_corrected, Time == 0 & source_plate == plate, select = c("GR24","Drug", "cell", "Final_conc_uM", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp <- tmp %>% 
    dplyr::group_by(Drug, cell,Final_conc_uM) %>% 
    dplyr::summarise(median_GR24 = median(GR24))
  
  lapply(unique(tmp$Drug), function(drug_idx){
    tmp <- subset(tmp, Drug == drug_idx)
    tmp$Final_conc_uM <- tmp %>%  dplyr::group_indices(Final_conc_uM)
    return(tmp)
  }) -> tmp
  
  tmp <- do.call(rbind,tmp)
  
  control <- subset(tmp, Drug == "DMSO")
  
  drug <- subset(tmp, Drug %in% drugs_in_screen)
  
  tmp <- dplyr::right_join(drug, control[,c("cell", "median_GR24")], by = "cell")
  
  tmp$percent_change_GR <- (tmp$median_GR24.x/tmp$median_GR24.y) *100
  
  return(tmp)
})-> output_GR24

output_GR24 <- do.call(rbind,output_GR24)

# save output 

setwd(path_data_file)

write.csv(output_GR24, "outcomes_GR24.csv")

# #calculate percent inhibition for confluence ----------------------------------

lapply(c("P1", "P2"), function(plate){ 
  # calculate results by plate
  #plate = "P1"
  tmp <- subset(data_corrected, Time >=0 & source_plate == plate, select = c("Conf",'Time',"Drug", "cell", "Final_conc_uM", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp <- tmp %>% 
    dplyr::group_by(Drug, cell,Final_conc_uM,Time) %>% 
    dplyr::summarise(median_conf = median(Conf))
  
  tmp <- tmp%>% dplyr::group_by(cell,Drug,Time)%>% dplyr::arrange(Final_conc_uM) %>% dplyr::mutate(Final_conc_uM = sequence(n()))
  
  tmp$cell_time <- paste(tmp$cell, tmp$Time)
  
  control <- subset(tmp, Drug == "DMSO")
  
  drug <- subset(tmp, Drug %in% drugs_in_screen)
  
  tmp <- dplyr::right_join(drug, control[,c("cell_time", "median_conf")], by = "cell_time")
  
  tmp$percent_change_GR <- (tmp$median_conf.x/tmp$median_conf.y) *100
  
  return(tmp)
})-> output_conf

output_conf <- do.call(rbind,output_conf)

# save output 

setwd(path_data_file)

write.csv(output_conf, "outcomes_conf_change_to_control.csv")

# plot phenotypic results -------------------------------------------------

#plot results for GR24 as heatmap

tmp <- subset(data_corrected,Time == 24, select= c('Drug', 'cell', "GR24", "Final_conc_uM"))

tmp <- tmp%>% dplyr::group_by(cell, Drug, Final_conc_uM) %>% dplyr::summarise(GR24 = mean(GR24))

tmp <- tmp%>% group_by(cell, Drug) %>%
  dplyr::slice(which.max(Final_conc_uM))

tmp$Final_conc_uM <- NULL

tmp <- pivot_wider(tmp, values_from = GR24, names_from = "cell")

tmp_names <- tmp$Drug

tmp$Drug <- NULL

rownames(tmp) <- tmp_names

setwd(paste(path_fig,"GR24", sep = "\\"))
png(paste("GR24",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey', scale = 'col')
dev.off()

#plot results for IC50 as heatmap 

tmp <- output_IC50[,c("cell_line", "agent","IC50")]

tmp$IC50 <- log10(tmp$IC50)

tmp <- pivot_wider(tmp, values_from = IC50, names_from = "cell_line")

tmp <- data.frame(tmp)

tmp_names <- tmp$agent

tmp$agent <- NULL

rownames(tmp) <- tmp_names

tmp <- t(apply(tmp, 1,function(x) ifelse(is.infinite(x),max(x[is.finite(x)]*1,na.rm = T),x)))

tmp <- ifelse(is.infinite(tmp),NA,tmp)

tmp <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]

labels_heatmap <- cell_metadata[order(match(cell_metadata$cell_line_display_name, colnames(tmp))),]

library(RColorBrewer)

labels_heatmap$colors <- as.character(factor(labels_heatmap$lineage_1, labels = RColorBrewer::brewer.pal(length(unique(labels_heatmap$lineage_1)), "Spectral")))

setwd(paste0(path_fig,'\\phenotype_figures'))

drugs_of_heatmap <- rownames(tmp)

pdf("IC50_phenotypes_72_legend.pdf")
heatmap_dendogram <- heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "RdBu"), margins=c(16,16),
          cexRow = 0.85,cexCol = 0.85,scale = 'row',ColSideColors = labels_heatmap$colors)

legend(x="bottomright", legend=c(labels_heatmap[!duplicated(labels_heatmap$lineage_1),]$lineage_1), 
       fill=labels_heatmap[!duplicated(labels_heatmap$lineage_1),]$colors)
dev.off()

#import GR24 data

setwd(path_data_file)

output_GR24 <- read.csv("outcomes_GR24.csv")

tmp <- output_GR24

col_breaks <- seq(0,120,length.out = 10)

col_idxs <- (RColorBrewer::brewer.pal(n=9, "RdYlBu"))

lapply(unique(tmp$Drug), function(drug_idx){
  #plot results by drug
  drug_data <- subset(tmp, Drug == drug_idx, select=c(cell, Final_conc_uM, percent_change_GR))
  
  threshold_low <- ifelse(drug_data$percent_change_GR<0,"<0","")
  
  threshold_high <- ifelse(drug_data$percent_change_GR>120,">120","")
  text_data <- cbind(drug_data[,c("cell", "Final_conc_uM")], data.frame(threshold_high))
  
  text_data$threshold_high <- paste0(text_data$threshold_high, threshold_low)
  
  drug_data$percent_change_GR <- ifelse(drug_data$percent_change_GR<0,0,drug_data$percent_change_GR)
  
  drug_data <- pivot_wider(drug_data, values_from = percent_change_GR, names_from = Final_conc_uM)
  colnames(drug_data) <- paste("conc=",colnames(drug_data))
  
  tmp_name <- drug_data$`conc= cell`
  
  drug_data$`conc= cell` <- NULL
  
  rownames(drug_data) <- tmp_name
  
  
  text_data <- pivot_wider(text_data, values_from = threshold_high, names_from = Final_conc_uM)
  colnames(text_data) <- paste("conc=",colnames(text_data))
  
  tmp_name <- text_data$`conc= cell`
  
  text_data$`conc= cell` <- NULL
  
  rownames(text_data) <- tmp_name
  
  setwd(paste(path_fig,'RSgroups', sep = "\\"))
  png(paste("drug=",drug_idx,"percent_change_GR24.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3)
  dev.off()
  
  png(paste("drug=",drug_idx,"percent_change_GR24_cluster.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3,Colv='none')
  dev.off()
  
  
})

# plot variance for each pair of cell and drug

GR24_var <- output_GR24 %>% dplyr::group_by(Drug,cell) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

GR24_var <- pivot_wider(GR24_var,values_from = variation, names_from = Drug)
GR24_var <- data.frame(GR24_var)
rownames(GR24_var) <- GR24_var$cell
GR24_var$cell <- NULL

labels_heatmap <- cell_metadata[order(match(cell_metadata$cell_line_display_name, rownames(GR24_var))),]


labels_heatmap$colors <- as.character(factor(labels_heatmap$lineage_1, labels = RColorBrewer::brewer.pal(length(unique(labels_heatmap$lineage_1)), "Spectral")))
library(RColorBrewer)

m=as.matrix((GR24_var))
plt <- heatmap(m,scale = "row",col = bluered(100), RowSideColors = labels_heatmap$colors)
dev.off()

# plot basal growth rates at time of sampling, drug == PBS

GR_basal_24h <- subset(data_corrected, Time == 24 & Drug %in% c("PBS"))

GR_basal_24h$cell <- factor(GR_basal_24h$cell, levels = labels_heatmap$cell_line_display_name[rev(plt$rowInd)])

ggplot(GR_basal_24h, aes(cell, GR24))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0,0.07))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

# plot avg basal growth rates at time of sampling, drug == PBS

GR24_basal_avg <- data.frame(GR_basal_24h %>% group_by(cell) %>% dplyr::summarize(avg_basal_GR24 = median(GR24)))

rownames(GR24_basal_avg) <- GR24_basal_avg$cell

GR24_basal_avg$cell <- NULL

GR24_basal_avg$tmp <- GR24_basal_avg$avg_basal_GR24


GR24_basal_avg <- GR24_basal_avg[order(rownames(GR24_basal_avg)),] 

GR24_basal_avg <- GR24_basal_avg[rownames(GR24_basal_avg)[heatmap_dendogram$colInd],]


setwd(paste0(path_fig,'\\phenotype_figures'))

pdf("GR_median_24h.pdf")
heatmap.2(as.matrix(GR24_basal_avg), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PiYG"), margins=c(16,16),
          cexRow = 0.85,cexCol = 0.85,scale = 'none',Rowv = "null",Colv = "null")

dev.off()

# plot basal growth rates at time of sampling, drug == all controls

GR_basal_24h <- subset(data_corrected, Time == 24 & Drug %in% c("PBS", "DMSO"))

GR_basal_24h$cell <- factor(GR_basal_24h$cell, levels = labels_heatmap$cell_line_display_name[rev(plt$rowInd)])

ggplot(GR_basal_24h, aes(cell, GR24, colour = Drug))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0,0.07))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))


# plot variance across concentratin for each pair of cell and drug

GR24_var <- output_GR24 %>% dplyr::group_by(Drug,Final_conc_uM) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

GR24_var <- pivot_wider(GR24_var,values_from = variation, names_from = Final_conc_uM)
GR24_var <- data.frame(GR24_var)
rownames(GR24_var) <- GR24_var$Drug
GR24_var$Drug <- NULL

GR24_var <-GR24_var[rev(na.omit(drugs_of_heatmap[heatmap_dendogram$rowInd],match(rownames(GR24_var) ))),] 

library(viridis)

setwd(paste0(path_fig,'\\phenotype_figures'))

pdf("variance_across_concentration.pdf")
heatmap.2(abs(as.matrix(GR24_var)), trace="none", key=T,col = viridis(20), margins=c(16,16),
          cexRow = 0.65,cexCol = 0.65,scale = 'row',Rowv = "null",Colv = "null")
dev.off()


# compute concentrations variability within drugs based on GR24 ----------------------------

#The goal for this is to:
  # Remove GR24 outliers
  # remove pairs of drug_conc that are inefective across all CCLs or that are too strong
  #low/high = drug_conc without/with effect across all CCLs

#import GR24 data

setwd(path_data_file)

output_GR24 <- read.csv("outcomes_GR24.csv")

#calculate CV
GR24_outliers_low <- output_GR24 %>% dplyr::group_by(Drug, Final_conc_uM) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

#label low
GR24_outliers_low$outliers <- ifelse(GR24_outliers_low$variation<=  cutoff_GR_conc_min_CV,"low", NA)

#save output low

setwd(path_data_file)

write.csv(GR24_outliers_low, "GR24_outliers_low.csv")

#define high outliers

GR24_outliers_high <- output_GR24

GR24_outliers_high$outliers <- ifelse(GR24_outliers_high$percent_change_GR<=  cutoff_GR_max_growth_effect,"high", NA)

#save output

setwd(path_data_file)

write.csv(GR24_outliers_high, "GR24_outliers_high.csv")


# plot exclusions  

exclusions_min_growth <-GR24_outliers_low[GR24_outliers_low$variation<=  cutoff_GR_conc_min_CV,]

tmp <- exclusions_min_growth

tmp$outliers <- NULL

tmp$X <- NULL

tmp$variation <- log10(abs(tmp$variation))

tmp <- pivot_wider(tmp, values_from = variation, names_from = Final_conc_uM)

colnames(tmp) <- paste("conc=",colnames(tmp))

tmp_name <- tmp$`conc= Drug`

tmp$`conc= Drug` <- NULL

rownames(tmp) <- tmp_name

exclusions <- ifelse(tmp <= log10(cutoff_GR_conc_min_CV), "off","")

setwd(paste(path_fig,'RSgroups', sep = "\\"))
png(paste("percent_change_CV.png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = rev(RColorBrewer::brewer.pal(n=11, "RdYlBu")), margins=c(16,16),
          cexRow = 1.5,cexCol = 2.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',cellnote = exclusions,notecol="black",notecex=2)
dev.off()