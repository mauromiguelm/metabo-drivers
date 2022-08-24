## This code analyse metabolomics data

# load packages and definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data_mean'
#path_data_file = "C:\\Users\\mauro\\Documents\\phd_results\\log2fc_full"
#path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures_mean\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\plate_converter.R")
library(openxlsx)
library(dplyr)
library(rhdf5)
library(ggplot2)
library(parallel)
library(proxy)
library(KEGGREST)
library(viridis)

substrRight <- function(x, n){
  
  if(nchar(x)>n){
    
    #x = substr(x, nchar(x)-n+1, nchar(x))
    
    #paste0("/",x)
    
    return("")
    
  }else{
    return(x)
    
  }
  
}


#import drug sensitity metrics

setwd(path_data_file)

data_GR50 <- read.csv("outcomes_growth_inhibition50.csv")
GR24_outliers_high <- read.csv('GR24_outliers_high.csv')
GR24_outliers_low <- read.csv('GR24_outliers_low.csv')

#import metabolomics

setwd(path_metabolomics_in)

dataContent<- h5ls("MM4_Mean mean_norm_DATA.h5")

data<- rhdf5::h5read(file = "MM4_Mean mean_norm_DATA.h5", '/data')

ions <- rhdf5::h5read(file = "MM4_Mean mean_norm_DATA.h5", '/annotation')

ions <- data.frame(ions)

#import cleaned metadata

setwd(paste(path_data_file,'metabolomics', sep = "//"))

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

# export log2fc data for clustering -------------------------------------

#map concentrations

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

lapply(unique(metab_fcs$drug), function(drug_idx){
  tmp <- subset(metab_fcs, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> metab_fcs

metab_fcs <- do.call(rbind, metab_fcs)

#remove concs wih no drug effect or too strong drug effect

metab_fcs$cell_drug_conc <- paste(metab_fcs$cell_line,metab_fcs$drug,metab_fcs$conc)

setwd(path_data_file)

groups_to_keep <-read.csv('outcomes_GR24.csv')
#
# groups_to_keep <- subset(RS_groups, !is.na(percent_change_GR))
#
groups_to_keep$cell_drug_conc <- paste(groups_to_keep$cell,groups_to_keep$Drug,groups_to_keep$Final_conc_uM)

#metab_fcs <- subset(metab_fcs, cell_drug_conc %in% unique(groups_to_keep$cell_drug_conc))

#return a matrix of data (rows are samples, cols are ion features) and metadata (sample, growth inhibition)

setwd(path_data_file)

lapply(unique(metab_fcs$drug), function(drug_idx){
  #drug_idx = "Methotrexate"
  print(drug_idx)
  tmp <- subset(metab_fcs, drug == drug_idx)

  tmp$pvalue <- NULL

  tmp <- tidyr::spread(tmp, key = ionIndex,value = log2fc)

  metadata <- tmp[,1:6]

  data <- tmp[,7:ncol(tmp)]

  metadata <- merge(metadata, groups_to_keep[,c('cell_drug_conc','percent_change_GR')], by ='cell_drug_conc',all.x = T)

  write.csv(metadata, paste("metadata",drug_idx,"log2fc.csv",sep ="_"))
  write.csv(data, paste("data",drug_idx,"log2fc.csv",sep ="_"))
})



# export log2fc data for tdsR -------------------------------------

#map concentrations

setwd(paste0(path_data_file,"\\tdsr"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

lapply(unique(metab_fcs$drug), function(drug_idx){
  tmp <- subset(metab_fcs, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> metab_fcs

metab_fcs <- do.call(rbind, metab_fcs)

#remove concs wih no drug effect or too strong drug effect

ions_of_interest <- which(ions$name %in%c("Asparagine", "Aspartate"))

ions_of_interest <- ions[ions_of_interest,]

metab_tds <- subset(metab_fcs, ionIndex %in% ions_of_interest$ionIndex & drug == "Asparaginase")

write.csv(metab_tds, "metabolomics-tdsr.csv")

# metabolome similarity across drugs_conc_cell lines ----------------------
#TODO replace raw intensities with log2fc when calculating simil

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

metab_fcs <- metab_fcs[metab_fcs$drug %in% unique(metab_fcs$drug)[1:10],]

groups <- unique(paste(metab_fcs$cell, metab_fcs$drug, metab_fcs$conc,sep = '_'))

groups <- t(combn(groups,m = 2))

metab_fcs$groups <- paste(metab_fcs$cell, metab_fcs$drug, metab_fcs$conc,sep = '_')

numWorkers <- 7

cl <-makeCluster(numWorkers, type="PSOCK")

calculate_slope_basal <- function(idx, metadata_, groups_) {
  #idx = 1
  meta_g1 <- metadata_[metadata_$groups %in% groups_[idx,1],]
  meta_g2 <- metadata_[metadata_$groups %in% groups_[idx,2],]

  data_g1 <- meta_g1$log2fc
  data_g2 <- meta_g2$log2fc

  #TODO release a warning in case data is empty

  #average groups

  #measures of similarity cosine

  d_cos <- proxy::simil(rbind(data_g1,data_g2),method = 'cosine')[1]
  #measures of similarity pearson

  d_pearson <- cor(data_g1,data_g2,method = 'pearson')
  #measures of similarity euclidean

  d_euc <- proxy::simil(rbind(data_g1,data_g2),method = 'euclidean')[1]
  #measures of similarity Hamman

  d_ham <- proxy::simil(rbind(data_g1,data_g2),method = 'Hamman')[1]

  #measures of similarity spearman

  d_spe <- cor(data_g1,data_g2,method = 'spearman')

  return(data.frame(g1 = groups_[idx,1], g2 = groups_[idx,2],
                    cosine = d_cos,
                    pearson = d_pearson,
                    euclidean = d_euc,
                    hamman = d_ham,
                    spearman = d_spe))

}

system.time({
  results <- parallel::parLapply(cl=cl,1:nrow(groups), calculate_slope_basal, metadata_=metab_fcs, groups_ = groups)
})

parallel::stopCluster(cl)

rm(cl)
# save results

results <- do.call(rbind, results)

setwd(paste0(path_data_file,"\\metabolomics","\\similarity"))

write.csv(results, 'similarity_results_log2fc.csv')

# plot results

setwd(paste0(path_data_file,"\\metabolomics","\\similarity"))

results <- read.csv('similarity_results_log2fc.csv')

results <- subset(results,grepl(results[,2],pattern = 'Erlotinib'))

results <- subset(results,grepl(results[,3],pattern = 'Erlotinib'))

tmp <- results

tmp2 <- results

tmp2$g2 <- tmp$g1
tmp2$g1 <- tmp$g2

tmp <- rbind(tmp, tmp2)

wide_result <- tidyr::spread(tmp[,c('g1','g2','cosine')],'g2','cosine',fill = 1)

rownames(wide_result) <- wide_result$g1

wide_result$g1 <- NULL
col <- viridis::viridis(999)

setwd(path_fig)
library(heatmaply)

png('heatmap_cos_sim.png',width = 1000,height = 1000)

heatmaply(as.matrix(wide_result))
plt <- gplots::heatmap.2(as.matrix(wide_result),symm = T,trace = 'none', col=col)
dev.off()

# plot one drug

# make a correlation network

require(xts)
require(quantmod)
require(igraph)
# cor_mat<- matrix( runif(100), nr=10 )

cor_mat <- wide_result

cor_mat[ lower.tri(cor_mat, diag=TRUE) ]<- 0

cor_mat[ abs(cor_mat) < 0.6]<- 0

graph <- graph.adjacency(abs(cor_mat)>0.6, mode="upper",weighted = T)

#E(graph)$weight<-t(cor_mat)[abs(t(cor_mat))>0.7]

E(graph)$weight<-t(cor_mat)[abs(t(cor_mat))>0.6]

v_names <- sapply(strsplit(V(graph)$name,split = '_'),"[[",2)

discrete_colors <- viridis::viridis(length(unique(v_names)))

v_colors <- factor(v_names,labels = discrete_colors)

V(graph)$color <- v_colors

graph$layout <- layout.circle

map <- data.frame(v_names,v_colors)

map <- map[!duplicated(map[,1]),]

graph$layout <- layout.fruchterman.reingold

png('corr_network.png',width=500)
plot(decompose.graph(graph)[[which.max(sapply(decompose.graph(graph), vcount))]],frame=T,
     edge.arrow.size=0.5,
     vertex.label.cex=0.7,vertex.size=3, col = V(graph)$color)
legend('topleft',title="Colors", cex=0.75, pch=16,
       col=map[,2],
       legend=map[,1], ncol=2)
dev.off()


# association GR24 with log2fc --------
# lm(metab~GR50) across all cell lines and concentrations that we have filtered strong effects/unnefective concentrations
# import log2fc data

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

lapply(unique(metab_fcs$drug), function(drug_idx){
  tmp <- subset(metab_fcs, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> metab_fcs

metab_fcs <- do.call(rbind, metab_fcs)

#for every drug, run 100.000 iterations
#distribute these iterations across all ions

bootstrap_association_ion_drug_effect <- function(drug_idx,fc_data,pheno_data, nboot=100000){
  #drug_idx = 'Oxfenicine'

  #distribute nboot across all ions
  print(paste('drug:',paste(drug_idx),'is not running'))
  tmp <- subset(fc_data, drug == drug_idx)

  tmp$cell_conc <- paste(tmp$cell_line, tmp$conc,sep = "_")

  tmp_growth_metrics <- subset(pheno_data, Drug == drug_idx, select=c("cell",'Final_conc_uM', 'percent_change_GR'))

  tmp_growth_metrics <- subset(tmp_growth_metrics, !is.na(percent_change_GR))

  tmp_growth_metrics$cell_conc <- paste(tmp_growth_metrics$cell, tmp_growth_metrics$Final_conc_uM,sep = "_")

  tmp <- merge(tmp, tmp_growth_metrics, by = "cell_conc")

  nions <- length(unique(tmp$ionIndex))

  iter_nr <- round(nboot/nions)

  iter_per_ion <- rep(iter_nr,nions)

  residual <- sum(iter_per_ion)-nboot


  if(length(unique(tmp$cell_conc))>10){

    out_boot <- data.frame()
    #make sure df is not empty

    #correct iterations for the residual of nboot

    iter_per_ion[sample(1:nions, abs(residual))] <- iter_nr - residual/abs(residual)

    stopifnot(sum(iter_per_ion)==nboot)

    #create permutations

    perm_df <- t(replicate(iter_nr, sample(length(unique(tmp$cell_conc)), replace = FALSE, prob = NULL)))

    #for every ion in each drug associate with GR50



    for(ion_idx in 1:nions){
      #ion_idx=1
      print(paste(drug_idx,ion_idx))

      tmp_ion <- subset(tmp, ionIndex == ion_idx)

      count_iter <- 0

      while(count_iter < iter_per_ion[ion_idx]){

        slope <- lm(log2fc~log10(tmp_ion$percent_change_GR[perm_df[count_iter+1,]]), tmp_ion)

        pvalue <- summary(slope)

        out_boot <- rbind(out_boot, data.frame(ion_idx,slope$coefficients[[2]],pvalue$r.squared,pvalue$adj.r.squared,pvalue$coefficients[2,4]))

        count_iter = count_iter +1

        }
    }
    print(paste('drug:',paste(drug_idx),'job finished'))

    return(out_boot)
  }

}

numWorkers <- 7

cl <-makeCluster(numWorkers, type="PSOCK",outfile = "tmp_err_boot-association.txt")

system.time({
  results <- parallel::parLapply(cl=cl,unique(metab_fcs$drug), bootstrap_association_ion_drug_effect, fc_data=metab_fcs,
                                 pheno_data = RS_groups)
})

parallel::stopCluster(cl)

rm(cl)

#save results bootstrap

setwd(paste0(path_data_file,"\\metabolomics\\log2fc"))

names(results) <- unique(metab_fcs$drug)

save(results,file = 'bootstrap_results.Rdata')

#report bootstrap confidence interval as table

lapply(names(results), function(drug_idx){

  #drug_idx = 'Methotrexate'
  tmp <- results[[drug_idx]]


  #get 0.025 and 0.975 quantiles which correcponds to 95% CI for the distribution
  conf_intervals <- t(data.frame(quantile(tmp$slope.coefficients..2..,probs = c(0.005, 0.995))))

  rownames(conf_intervals) <- drug_idx

  #plot results

  #hist(tmp$slope.coefficients..2..,breaks = 1000)
  #abline(v=conf_intervals[1],col="blue",lwd=2)
  #abline(v=conf_intervals[2],col="red",lwd=2)

  return(conf_intervals)

})-> bootstrap_ci

bootstrap_ci <- do.call(rbind, bootstrap_ci)

#save confidence interval results

setwd(paste0(path_data_file,"\\metabolomics\\log2fc"))

write.csv(bootstrap_ci,file = 'confidence_intervals_from_permutation_1.csv')

#calculate metrics for each drug_ions of log2fc vs. GR24

lapply(unique(paste(metab_fcs$drug,metab_fcs$ionIndex, sep = "_")),function(drug_ion){
  #for every ion in each drug associate with GR50
  print(drug_ion)
  #drug_ion <-  "Decitabine"

  tmp <- subset(metab_fcs, drug == strsplit(drug_ion, split = "_")[[1]][1] &
                ionIndex == as.numeric(strsplit(drug_ion, split = "_")[[1]][2]))

  tmp$conc <- tmp %>% dplyr::group_by(concentration) %>%  dplyr::group_indices(concentration)
  tmp$cell_conc <- paste(tmp$cell_line, tmp$conc,sep = "_")


  tmp_growth_metrics <- subset(RS_groups, Drug == strsplit(drug_ion, split = "_")[[1]][1], select=c("cell",'Final_conc_uM', 'percent_change_GR'))

  tmp_growth_metrics <- subset(tmp_growth_metrics, !is.na(percent_change_GR))

  tmp_growth_metrics$cell_conc <- paste(tmp_growth_metrics$cell, tmp_growth_metrics$Final_conc_uM,sep = "_")

  if(nrow(tmp_growth_metrics) > 10){
    #only include drugs where n of defined G50 is higher or equal than 3

    #combine metabolomics with GR50

    tmp <- merge(tmp, tmp_growth_metrics, by = "cell_conc")

    #calculate linear regression intensity vs. GI50

    slope <- lm(log2fc~log10(percent_change_GR), tmp)

    pvalue <- summary(slope)

    return(c(drug_ion,slope$coefficients[2],pvalue$r.squared,pvalue$adj.r.squared,pvalue$coefficients[2,4]))


  }else{
    return(NA)
  }

}) -> slope_metabolite_effect_on_growth

slope_metabolite_effect_on_growth <- do.call(rbind, slope_metabolite_effect_on_growth)

slope_metabolite_effect_on_growth <-data.frame(slope_metabolite_effect_on_growth)

slope_metabolite_effect_on_growth$log10.percent_change_GR. <- as.numeric(slope_metabolite_effect_on_growth$log10.percent_change_GR.)
slope_metabolite_effect_on_growth$V4 <- as.numeric(slope_metabolite_effect_on_growth$V4)
slope_metabolite_effect_on_growth$V5 <- as.numeric(slope_metabolite_effect_on_growth$V5)


slope_metabolite_effect_on_growth <- cbind(stringr::str_split_fixed(string = slope_metabolite_effect_on_growth$V1,"_",n = 2),slope_metabolite_effect_on_growth[,-c(1)])

slope_metabolite_effect_on_growth <- data.frame(slope_metabolite_effect_on_growth)

colnames(slope_metabolite_effect_on_growth) <- c('drug', 'ionIndex', 'slope', 'r2', 'adj-r2', 'pvalue')

#save associations

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

write.csv(slope_metabolite_effect_on_growth, 'metabolite_GR24_association.csv')

#select the ions that survived permutation thresohld

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

df <- read.csv('metabolite_GR24_association.csv')

lapply(unique(df$drug), function(drug_idx){

  tmp <- subset(df, drug == drug_idx & !is.na(slope))

  if(nrow(tmp)>1){
    tmp_ci <- bootstrap_ci[drug_idx,]

    tmp <- subset(tmp,slope < tmp_ci[1] | slope> tmp_ci[2])

  }
  return(tmp)
}) -> df

df <- do.call(rbind,df)

#plot number of associations per drug

df%>% group_by(drug) %>% summarise(association_count = n()) %>%
  ggplot(aes(x=reorder(drug,association_count), y=association_count ))+
  geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab('Drug')+
  ylab('Number of associations')

df$X <- seq(1,nrow(df))

rownames(df) <- NULL

colnames(df)[1] <- "idx"

#save associations

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

write.csv(df, 'metabolite_GR24_association_survivingCI.csv')

# baseline vs drug treated ion compairison --------------------------------

lapply(unique(ions$ionIndex), function(ionidx){
  #ionidx = 816
  print(ionidx)
  #iterate over ions and claculate association betwen basal and drug
  lapply(unique(RS_groups$Drug), function(drugidx){
    #drugidx = 'Methotrexate'
    #iterate over drug

    RSgroup_sub <- subset(RS_groups, Drug == drugidx & !is.na(percent_change_GR))

    metadata_doi <- subset(metadata, drug == drugidx)

    metadata_doi$conc <- metadata_doi %>% dplyr::group_by(conc) %>%  dplyr::group_indices(conc)

    metadata_doi <- subset(metadata_doi, conc %in% unique(RSgroup_sub$Final_conc_uM))

    metadata_control <- subset(metadata, cell %in% unique(metadata_doi$cell) & source_plate == unique(metadata_doi$source_plate) & drug == 'DMSO' &
                                 conc == 367)

    if(nrow(metadata_doi)>7){
      #skip drugs with low number of R/S cell lines

      if(!length(unique(metadata_doi$cell))==length(unique(metadata_control$cell))){stop('diverging number of cell lines across groups')}

      data_median_control <- data[ionidx,metadata_control$idx]

      data_median_control <- cbind(metadata_control,data.frame("intensities" = data_median_control))
  
      data_median_control <- data_median_control %>% dplyr::group_by(cell) %>% dplyr::summarize(median_ion = median(intensities))

      data_median_drug <- data[ionidx,metadata_doi$idx]

      data_median_drug <- cbind(metadata_doi,data.frame("intensities" = data_median_drug))
      
      data_median_drug <- data_median_drug %>% dplyr::group_by(cell) %>% dplyr::summarize(median_ion = median(intensities))

      comb_data <- merge(data_median_control, data_median_drug, by = 'cell')
      
      slope <- lm(median_ion.y~median_ion.x, comb_data)

      pvalue <- summary(slope)

      return(c(drugidx,ionidx,slope$coefficients[2],pvalue$r.squared,pvalue$adj.r.squared,pvalue$coefficients[2,4]))

    }

  })


}) -> basal_drug_association

basal_associations <- unlist(basal_drug_association, recursive = F)

basal_associations <- do.call(rbind, basal_associations)

basal_associations <-data.frame(basal_associations)

colnames(basal_associations) <- c('drug', 'ionIndex', 'slope', 'r2', 'adj-r2', 'pvalue')

#save results from basal association

#save associations

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

write.csv(basal_associations, 'basal_association.csv')

# Determining half effect concentration for relevant ions -----------------

#import GR24 associations with metabolism

#read associations

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

df <- read.csv('metabolite_GR24_association_survivingCI.csv')

df$X <- NULL

#import full GR24 data

setwd(paste0(path_data_file))

outcomes_GR24 <- read.csv("outcomes_GR24.csv")

#for each ion that is associated with growth, use all concentrations
#to determine IC50

drug_cell_ED50 = list()
models_GR24 = list()
fixed_u = 1

scale_y = function(x){(x-min(x))/(max(x)-min(x))}

for(drug in unique(df$drug)){
  #print(drug)

  #drug = "YC-1"
  #drug = "Decitabine"

  sub_drug = subset(outcomes_GR24, Drug == drug)
  
  sub_drug$percent_change_GR <- ifelse(sub_drug$percent_change_GR>100,100, sub_drug$percent_change_GR )

  #calculate ED50 for every cell line
  sub_drug$norm_change_GR <- scale_y(sub_drug$percent_change_GR)

  for(cell in unique(sub_drug$cell)){
    print(cell)
    #only use drugs that have at least 40% reduction in GR to calculate
    #cell ='IGROV1'
    #cell = "BT549"
    sub_drug_cell <- sub_drug[sub_drug$cell == cell,]
    if(any(sub_drug_cell$percent_change_GR <2000)){

      #order to eliminate bias from LOQ, non-linear effects
      sub_drug_cell <- sub_drug_cell[order(sub_drug_cell$Final_conc_uM),]
      sub_drug_cell$norm_change_GR <- sort(sub_drug_cell$norm_change_GR,decreasing = T)

      #fit model
      model = try(drc::drm(norm_change_GR~Final_conc_uM,data = sub_drug_cell, fct=drc::l4(fixed = c(NA,NA, fixed_u, NA)),
                           na.action = na.omit), silent = T)

      #get params, ED50

      if(class(model)=='try-error'){
        h <-  NA
        E0 <- NA
        Einf <- NA
        ED50 <-  NA
        model = "undefined"
      }else{
        h <-  model$coefficients[1]
        E0 <- model$coefficients[2]
        Einf <- fixed_u
        ED50 <- model$coefficients[3]
  

      }


      drug_cell_ED50=append(drug_cell_ED50,list(data.frame(
        drug,
        cell,
        h,
        E0,
        Einf,
        ED50
        )))
      
      

      models_GR24 = append(models_GR24, list(model))

    }
  }
}

drug_cell_ED50 <- do.call(rbind, drug_cell_ED50)

drug_cell_fc <- list()
models_log2fc <- list()

for(idx in df$idx){
  print(idx)
  #idx = 492 #hbp
  #idx = 138 orotate
  sub_df <- df[idx,]

  sub_fcs <- subset(metab_fcs, drug == sub_df$drug & ionIndex==sub_df$ionIndex)
  
  sub_fcs$norm_log2fc <- scale_y(abs(sub_fcs$log2fc))

  for(cell in unique(sub_fcs$cell_line)){
    #cell = 'A498'
    sub_cell_fcs <- subset(sub_fcs, cell_line == cell)
    #order to eliminate bias from LOQ, non-linear effects
    sub_cell_fcs <- sub_cell_fcs[order(sub_cell_fcs$conc),]
    sub_cell_fcs$norm_log2fc <- sort(sub_cell_fcs$norm_log2fc,decreasing = T)

    #fit model
    
    model = try(drc::drm(norm_log2fc~conc,data = sub_cell_fcs, fct=drc::l4(fixed = c(NA,NA, fixed_u, NA)),
                         na.action = na.omit), silent = T)

    #get params, ED50
    
    if(class(model)=='try-error'){
      h <-  NA
      E0 <- NA
      Einf <- NA
      ED50 <- NA
      model <- "Undefined"
    }else{
      h <-  model$coefficients[1]
      E0 <- model$coefficients[2]
      Einf <- fixed_u
      ED50 <- model$coefficients[3]
      
      
    }
    
    drug_cell_fc=append(drug_cell_fc,list(data.frame(
      idx,
      cell,
      h,
      E0,
      Einf,
      ED50
    )))

    
    models_log2fc = append(models_log2fc, list(model))

  }
}


drug_cell_fc <- do.call(rbind, drug_cell_fc)

#save results

setwd(paste(path_data_file, "GI50_associations", sep ="\\"))

write.csv(drug_cell_fc, "drug_cell_fc.csv")

write.csv(drug_cell_ED50,"drug_cell_ED50.csv")

#TODO save models, so I can use them later

#read results

setwd(paste(path_data_file, "GI50_associations", sep ="\\"))

drug_cell_fc <- read.csv("drug_cell_fc.csv")

drug_cell_ED50 <- read.csv("drug_cell_ED50.csv")

#gap fill NA with high concentration

lapply(unique(drug_cell_fc$idx), function(idx){
  #idx = 10
  sub_data <- drug_cell_fc[drug_cell_fc$idx==idx,]
  
  sub_data$ED50 <- ifelse(is.na(sub_data$ED50), max(sub_data$ED50, na.rm = T)*1.1, sub_data$ED50)
  
  return(sub_data)
}) -> tmp


setwd(paste(path_data_file, "GI50_associations", sep ="\\"))

library(ggplot2)

lapply(tmp, function(sub_data){
  #sub_data <- tmp[[600]]
  
  #compare ED50 for log2fc and GR24
  sub_fc <- df[which(df$idx==unique(sub_data$idx)),]
  print(unique(sub_data$idx))
  sub_drug <- subset(drug_cell_ED50, drug == unique(sub_fc$drug))
  
  
  
  data_merged <- merge(sub_drug[,c("cell","ED50")], sub_data[,c("cell","ED50")], by = "cell")
  colnames(data_merged) <- c('cell',"GR24", "log2fc")
  data_merged <- tidyr::pivot_longer(data_merged, -cell)
  
  mean_data <- aggregate(value ~ name, data = data_merged, FUN= "median" )
  
  GR_diff = mean_data[mean_data$name=='GR24',2]/mean_data[mean_data$name=='log2fc',2]
  
  metab_name <- gsub(":", "-",sub_fc$name)
  metab_name <- gsub("/", "-",metab_name)
  
  if(abs(GR_diff)>1.5){
    ggplot(data_merged, aes(name, value, fill = name))+
      geom_boxplot()+
      # geom_point() is used to plot data points on boxplot
      geom_jitter(aes(fill=name,group=cell),size=3,shape=21)+
      ylab('ED50')+
      ggtitle(paste(sub_fc$drug, sub_fc$name, sep = "_"))-> plot
    ggsave(filename = paste(sub_fc$drug,metab_name,".png", sep ="_"), plot = plot)
    
  }
  
  return(data.frame(unique(sub_data$idx),GR_diff))
  
}) -> tmp

EC50_statistics <- data.frame(do.call(rbind,tmp))

colnames(EC50_statistics) <- c("association_index","EC50_FC")

setwd(paste(path_data_file, "GI50_associations", sep ="\\"))

write.csv(EC50_statistics, "statistics_EC50.csv")

library(ggplot2)

#plot fitting

plot(0,0,xlim = c(0,5),ylim = c(0,1),type = "n",xlab='Dose',ylab='Effect')
cols = rev(heat.colors(length(models)))
dose=seq(0,5,0.1)

for(x in 1:length(models_GR24)){
  print(x)

  model = models_GR24[[x]]

  if(class(model)=='drc'){
    tmp <- predict(model, data.frame(dose=dose),
                   interval = "prediction")
    lines(dose, tmp[,1],col= 'blue')

  }

  model = models_log2fc[[x]]

  if(class(model)=='drc'){
    tmp <- predict(model, data.frame(dose=dose),
                   interval = "prediction")
    lines(dose, tmp[,1],col= 'red')

  }
}



# combine all metrics into one matrix -------------------------------------


#read EC50 results

setwd(paste(path_data_file, "GI50_associations", sep ="\\"))

drug_cell_fc <- read.csv("drug_cell_fc.csv")

drug_cell_ED50 <- read.csv("drug_cell_ED50.csv")

EC50_statistics <- read.csv("statistics_EC50.csv")

EC50_statistics$X <- NULL

#select the ions that survived permutation thresohld

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

df <- read.csv('metabolite_GR24_association_survivingCI.csv')

df$X <- NULL

df$color <- ifelse(df$slope<0, "Negative", "Positive")

colnames(df)[4:8] <- paste("DMA", colnames(df)[4:8], sep = "_") 

tmp_ions <- ions[,c("ionIndex",'mzLabel','idKEGG','score', "name")] %>% dplyr::group_by(ionIndex) %>% dplyr::arrange(score) %>% dplyr::slice(n())

tmp_ions[tmp_ions$mzLabel == 'mz133.0145',"name"] <- "Malate"

df <- merge(df, tmp_ions, by='ionIndex')

write.csv(df, file = 'metabolite_GR24_association_survivingCI_ionLabel.csv')

# combine log2fc/GI50 association with basal/treated metabolome association

df$association_index = df$idx

df$drugion <- paste(df$drug, df$ionIndex)

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

basal_associations <- read.csv('basal_association.csv')

basal_associations$X <- NULL

basal_associations$drugion <- paste(basal_associations$drug, basal_associations$ionIndex)

colnames(basal_associations)[3:6] <- paste("basal", colnames(basal_associations)[3:6], sep = "_")

df <- merge(df, basal_associations[,c('drugion', 'basal_slope','basal_pvalue',"basal_r2", "basal_adj.r2")], by='drugion')

#combine GR50 association

colnames(EC50_statistics)[2] <- c("EC50_FC")

df <- merge(df, EC50_statistics, by= "association_index")

df$DMA <- ifelse(df$DMA_slope<0, "Sensitivity", "Resistant")
df$Basal_effect <- ifelse(df$basal_slope<1, "Consumption", "Synthesis")
df$GR_independence <- ifelse(df$EC50_FC >1, T, F)

#save full with all 3 metrics

setwd(path_data_file)
write.csv(df, "drug_metabolite_associations.csv")

#generate statistics as results per group

setwd(path_data_file)
df <- read.csv("drug_metabolite_associations.csv")


stats_associations <-  xtabs(~Basal_effect+DMA+GR_independence, df)

#generate number of associations per ion

stats_count_per_ion <- df %>% group_by(name) %>% summarize(count = n())


# plot associations by drug, heatmap by metrics ---------------------------

#import results

setwd(path_data_file)

df <- read.csv("drug_metabolite_associations.csv")
df$X = NULL

lapply(unique(df$drug), function(drug){
  #drug = "Decitabine"
  tmp <- df[df$drug==drug,]
  
  n_interactions <- length(unique(tmp$name))
  
  #subset ion and the tree metrics
  
  
  tmp <- tmp[,c("name","mzLabel","EC50_FC","DMA_slope","basal_slope")]
  
  tmp$name <- as.character(sapply(tmp$name, substrRight, 25))
  
  tmp$name <- paste(tmp$name, paste0("(",tmp$mzLabel,")"), sep = " ")
  
  tmp$mzLabel <- NULL
  
  colnames(tmp) <- c("Metabolite","log2(GR ind)", "DMA", "Basal effect")
  
  tmp$`log2(GR ind)` <- log2(tmp$`log2(GR ind)`)
  
  #limit GRI at 4
  
  tmp$`log2(GR ind)` <- ifelse(tmp$`log2(GR ind)`>4, 4, tmp$`log2(GR ind)`)
  
  tmp$`log2(GR ind)` <- ifelse(tmp$`log2(GR ind)`<c(-4), c(-4), tmp$`log2(GR ind)`)
  tmp <- tidyr::pivot_longer(tmp,-Metabolite, names_to = "metric")
  
  tmp$metric <- factor(tmp$metric,levels = c("log2(GR ind)", "Basal effect", "DMA"))
  
  GRI <- tmp[which(tmp$metric=="log2(GR ind)"),'value'][[1]]
  
  GRI_max <- max(abs(GRI))
  
  GRI_colmap = colorRamp2(c(-GRI_max, 0, GRI_max), c("darkgreen", "grey", viridis(100)[100]))
  
  GRI_col <- GRI_colmap(GRI)
  
  tmp$color <- NA
  tmp$color <- ifelse(tmp$metric=="DMA" & tmp$value<0, "#FF0000", tmp$color)
  tmp$color <- ifelse(tmp$metric=="DMA" & tmp$value>0, "green4", tmp$color)
  tmp$color <- ifelse(tmp$metric=="Basal effect" & tmp$value<1, "blue", tmp$color)
  tmp$color <- ifelse(tmp$metric=="Basal effect" & tmp$value>1, "red", tmp$color)
  tmp[which(tmp$metric=="log2(GR ind)"),'color'] <- GRI_colmap(tmp[which(tmp$metric=="log2(GR ind)"),'value'][[1]])
  
  tmp$value <- ifelse(tmp$metric=="log2(GR ind)", 2,tmp$value)
  
  setwd(paste0(path_fig, "\\association_heatmap_by_drug"))
  ggplot(tmp, aes(metric, Metabolite))+
    geom_point(aes(size = abs(value), fill = color), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0,round(max(abs(tmp$value)))+1),breaks = round(seq(0,round(max(abs(tmp$value)))+1,length.out=5))) + 
    labs( x= "", y = "", size = "Slope", fill = "")  + 
    theme(legend.key=element_blank(), 
          axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
          axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
          legend.text = element_text(size = 10, face ="bold", colour ="black"), 
          legend.title = element_text(size = 12, face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.position = "right")+
    scale_fill_identity(aes(color)) -> plt
  
  ggsave(paste(drug, "h1.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.34,width = 7, limitsize = F)
  ggsave(paste(drug, "h0.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.4,width = 7, limitsize = F)
  ggsave(paste(drug, "h2.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.5,width = 7, limitsize = F)
  ggsave(paste(drug, "h3.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.3,width = 9, limitsize = F)
  
  #create a legend for GR independence
  
  GRI <- seq(c(-GRI_max), GRI_max, length.out = 7)
  
  GRI_col <- GRI_colmap(GRI)

  GRI <- round(GRI,0)
  
  pdf(paste(drug, "legend.pdf")) 
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x = 0, y = 1, fill = GRI_col, legend = GRI, title = "log2(GR independence)")
  dev.off()

  })


# plot top ranked interesting associations -----------------------

#number of associations per metabolite with diferent drugs (plot)

tmp <- df %>% group_by(ionIndex) %>% dplyr::summarise(count_metab = n())

setwd(path_fig)

ggplot(tmp, aes(x=count_metab, fill = factor(count_metab)))+
  geom_point(stat ="count", size = 4)+
  scale_fill_manual("legend", values = viridis::rocket(15))+
  theme_bw()+
  ylab("Number of metabolites")+
  xlab("Number of drugs")+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),
        axis.ticks = element_line(size = 1),axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20))+
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,4))+
  theme(legend.position = "none") -> plt

ggsave("number_metabs_associated_with_number_drugs.pdf",plt,device = "pdf",width=5, height=4,scale = 0.9)

# #sankey plot of positive associations -----------------------------------

#first col is drug, second is metab_pathway or any other meta for metabolites

library(networkD3)
library(dplyr)
library(tidyr)

#prep data for individual drug metab assoation

df_sankey <- df[,c('drug','name','DMA_slope')]

names(df_sankey) <- c("name", 'year1','value')

df_sankey$groups <- paste(df_sankey$name, df_sankey$year1)

df_sankey <- df_sankey[!df_sankey$year1 %in% df_sankey$name,] #remove metabolites with drug names

#create links and nodes

links <-
  df_sankey[,c(1,2)] %>%mutate(row = row_number()) %>%  # add a row id
  pivot_longer(-row, names_to = "column", values_to = "source") %>%  # gather all columns
  mutate(column = match(column, names(df))) %>%  # convert col names to col ids
  group_by(row) %>%
  mutate(target = lead(source, order_by = column)) %>%  # get target from following node in row
  ungroup() %>%
  filter(!is.na(target))  # remove links from last column in original data

links$groups <- paste(links$source, links$target)

links <- merge(links, df_sankey[,c('groups',"value")],by="groups")

nodes <- data.frame(name = unique(c(links$source, links$target)))
nodes$label <- sub('_[0-9]*$', '', nodes$name) # remove column id from node label

links$source_id <- match(links$source, nodes$name) - 1
links$target_id <- match(links$target, nodes$name) - 1
#links$value <- abs(links$value)

#plot

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source_id',
              Target = 'target_id', Value = 'value', NodeID = 'label',height = 4000,width = 1000)


# #sankey plot on pathways ------------------------------------------------

#prep data for individual drug to pathway association

#expand kegg id

df_sankey <- df[,c('drug','idKEGG')]

df_sankey <- tidyr::separate_rows(df_sankey,idKEGG, convert = TRUE, sep = ' ')

df_sankey <- subset(df_sankey,!grepl("^\\s*$", idKEGG))

# relate kegg compound ID to pathway

map_pathway_to_cpds <- data.frame(do.call(rbind,strsplit(temp, "\t")))

map_pathway_to_cpds$X1 <- gsub("path:","",map_pathway_to_cpds$X1)

hsa_metab_pathways <- data.frame(path_names = keggGet('path:hsa01100')[[1]][['REL_PATHWAY']])

which_pathways <- sapply(map_pathway_to_cpds, "%in%",rownames(hsa_metab_pathways))

map_pathway_to_cpds <- map_pathway_to_cpds[which_pathways,]

path_names <- lapply(df_sankey$idKEGG,function(idx){
  #idx =  "C00112"
  which_path <- apply(map_pathway_to_cpds,1,"%in%",idx)

  which_path <- apply(which_path, 1, any)

  if(sum(which_path)>0){


    path_names <- hsa_metab_pathways[rownames(hsa_metab_pathways)%in%map_pathway_to_cpds[which_path,1],]
    return(paste(path_names, collapse = "_"))
    }
  })

df_sankey$path_names <- path_names

df_sankey <- tidyr::separate_rows(df_sankey,path_names, convert = TRUE, sep = '_')

df_sankey <- df_sankey[,c('drug','path_names')]

names(df_sankey) <- c("name", 'year1')

library(tidyr)
library(dplyr)
#create links and nodes

links <-
  df_sankey %>%mutate(row = row_number()) %>%  # add a row id
  pivot_longer(-row, names_to = "column", values_to = "source") %>%  # gather all columns
  mutate(column = match(column, names(df))) %>%  # convert col names to col ids
  group_by(row) %>%
  mutate(target = lead(source, order_by = column)) %>%  # get target from following node in row
  ungroup() %>%
  filter(!is.na(target))  # remove links from last column in original data

nodes <- data.frame(name = unique(c(links$source, links$target)))
nodes$label <- sub('_[0-9]*$', '', nodes$name) # remove column id from node label

links$source_id <- match(links$source, nodes$name) - 1
links$target_id <- match(links$target, nodes$name) - 1
links$value <- 1

library(networkD3)

#plot

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source_id',
              Target = 'target_id', Value = 'value', NodeID = 'label')



# #single ion count per pathay --------------------------------------------

tmp <- df_sankey %>% group_by(name, year1) %>% summarise(count = n())

tmp <- subset(tmp, count >=4)

library(ggplot2)
ggplot(tmp, aes(y= name, x= year1, fill = count))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45,hjust =1))



#clirclize plot of top 10 associations ---------------------------------

setwd(path_data_file)
df <- read.csv("drug_metabolite_associations.csv")

df <- df %>% group_by(drug) %>% dplyr::arrange(abs(DMA_slope)) %>% dplyr::slice_tail(n=10)

df <- df %>% group_by(drug) %>%dplyr::arrange(abs(DMA_slope)) %>%  dplyr::mutate(slope_norm = c(10:1)[1:n()])

df$basal_slope <- ifelse(df$basal_slope >=2,2, df$basal_slope)
df$basal_slope <- ifelse(df$basal_slope <=c(0),-0, df$basal_slope)

df$DMA_slope <- ifelse(df$DMA_slope >=2,2, df$DMA_slope)
df$DMA_slope <- ifelse(df$DMA_slope <=c(-2),-2, df$DMA_slope)


df$name <- as.character(sapply(df$name, substrRight, 25))

df$name <- paste(df$name, paste0("(",df$mzLabel,")"), sep = " ")

df$EC50_FC <- log2(df$EC50_FC)

df$EC50_FC <- ifelse(df$EC50_FC >4,4, df$EC50_FC)
df$EC50_FC <- ifelse(df$EC50_FC <c(-4),-4, df$EC50_FC)

# merge MoA 

setwd(path_data_file)

drug_metadata <- read.csv("drug_metadata_moa.csv")

colnames(drug_metadata)[5] <- 'drug'

df <- merge(df, drug_metadata[,c('drug', 'Class')], by = 'drug')

#start circos plot



library(circlize)

df$x=0
df$color <- ifelse(df$DMA_slope<0, "#FF0000", "green4")
#circos.par("track.height" = 0.3,cell.padding = c(0.02, 0.04, 0.02, 0.04))
setwd(path_fig)

png("association_GR24_metabolite.png",width = 8000,height = 8000,res = 700)

x1 <- x2 <- y1 <-y2 <-0.5

circos.par("track.height" = 0.05,canvas.xlim = c(-(1+x1), 1+x2), canvas.ylim = c(-(1+y1), 1+y2),
           start.degree = 93)
circos.initialize(df$drug, x = (df$slope_norm))

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)

circos.trackText(sectors = df$drug,x = df$slope_norm, y= df$x, labels = df$name,
                 cex= 0.6,track.index = 1,facing = 'clockwise',col = df$color,adj =  c(0),niceFacing = T)

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)

circos.trackPoints(sectors = df$drug,x = df$slope_norm, y= df$x,cex = abs(df$DMA_slope)*0.5, col = df$color,pch = 16)

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)

col2 = colorRamp2(c(min(df$basal_slope), 1, max(df$basal_slope)), c("blue", "grey", "red"))


df$color2 <- col2(df$basal_slope)
circos.trackPoints(sectors = df$drug,x = df$slope_norm, y= df$x,col = df$color2 ,pch = 16)


# add EC50 analysis 

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)


GRI_max <- max(abs(df$EC50_FC))

col3 = colorRamp2(c(-GRI_max, 0, GRI_max), c("darkgreen", "grey", viridis(100)[100]))

df$color3 <- col3(df$EC50_FC)

circos.trackPoints(sectors = df$drug,x = df$slope_norm, y= df$x,col = df$color3 ,pch = 16)

# add drug MoA

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)

color4 <- RColorBrewer::brewer.pal(n = 9, name = 'Paired')
moa <- unique(df$Class)

df$color4 <- color4[match(df$Class, moa)]

i=1
for(x in (df$drug)){
  highlight.sector(track.index = 5,col = df$color4[i],sector.index = x)
  i=i+1
}

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,track.index = 6,
                           CELL_META$cell.ylim[2] - mm_y(13),
                           CELL_META$sector.index,cex = 1,facing = 'clockwise',niceFacing = T)
             },bg.border = NA)

dev.off()
circos.clear()



col3 = colorRamp2(c(-GRI_max, 0, GRI_max), c("darkgreen", "grey", viridis(100)[100]))

df$color3 <- col3(df$EC50_FC)


pdf(paste("circle-legend-GRI.pdf")) 
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = 0, y = 1, fill = col3(seq(-GRI_max,GRI_max,1)),legend = seq(-GRI_max,GRI_max,1), title = "log2(GR independence)")
dev.off()




pdf(paste("circle-legend-basal.pdf")) 
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = 0, y = 1, fill = col2(seq(0,2,0.25)),legend = seq(0,2,0.25), title = "Basal effect")
dev.off()



# plot summary of ED50


df$PASS <-  ifelse(df$EC50_FC>0,T,F)

EC50_stats_by_drug <- do.call(data.frame, aggregate(PASS~drug, df,  FUN = function(x) c(count = sum(x), total = length(x))))

EC50_stats_by_drug$PASS.percentage <- EC50_stats_by_drug$PASS.count/EC50_stats_by_drug$PASS.total
ggplot(EC50_stats_by_drug, aes(x=reorder(drug, PASS.percentage),y=PASS.percentage))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) -> plt

ggsave(paste0(path_fig,"","\\GI50_associations\\", "fraction_GI50.pdf"),plot = plt,  width = 25, height = 15, units = "cm")
  

# plot one interesting association ---------------------------------------------
#doi = drug of interest
#iot = ion of interest

doi <- "Decitabine"
iot <- 394

tmp2 <- subset(metab_fcs, ionIndex == iot & drug == doi)

lapply(unique(tmp2$drug), function(drug_idx){
  tmp <- subset(tmp2, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> tmp2

tmp2 <- do.call(rbind, tmp2)

tmp2$experiment <- paste(tmp2$cell_line, tmp2$drug,tmp2$conc, sep = "_")

tmp_growth_metrics <- subset(RS_groups, select=c('Drug',"cell",'Final_conc_uM', 'percent_change_GR'))

tmp_growth_metrics <- subset(tmp_growth_metrics, !is.na(percent_change_GR))

tmp_growth_metrics$experiment <- paste(tmp_growth_metrics$cell,tmp_growth_metrics$Drug, tmp_growth_metrics$Final_conc_uM,sep = "_")

tmp2 <- merge(tmp2, tmp_growth_metrics[,c('experiment', "percent_change_GR")], by = "experiment")

(ggplot(tmp2, aes(y=log2fc,x=log10(percent_change_GR), label = cell_line))+
  geom_point(col = 'red')+
  ggrepel::geom_text_repel(cex = 3)+
  theme_bw()+
  geom_smooth(method='lm', formula= y~x) -> plt)

plt <- plt+theme(axis.title =element_text(size=15))

setwd(paste(path_fig,'specific_examples', sep = "\\"))

ggsave(filename = paste(doi,iot,'.png', sep = "_"), plt,width = 4,height = 3.5,dpi = 300)

# dependency on concentration
#do increasing concentration change levels of these metabolites

ggplot(tmp2, aes(x = factor(conc), y=cell_line, fill = log2fc))+
  geom_tile()+
  scale_fill_gradient2(low = "red",mid = 'white', high = "blue",midpoint = 0)+
  theme_bw()+
  theme(axis.title =element_text(size=15), axis.text = element_text(size = 15))

# plot associations wiht basal metabolism

doi <- "Decitabine"
iot <- 394

RSgroup_sub <- subset(RS_groups, Drug == doi & !is.na(percent_change_GR))

metadata_doi <- subset(metadata, drug == doi)

metadata_doi$conc <- metadata_doi %>% dplyr::group_by(conc) %>%  dplyr::group_indices(conc)

metadata_doi <- subset(metadata_doi, conc %in% unique(RSgroup_sub$Final_conc_uM))

metadata_control <- subset(metadata, cell %in% unique(metadata_doi$cell) & source_plate == unique(metadata_doi$source_plate) & drug == 'DMSO' &
                             conc == 367)

data_mean_control <- data[iot,metadata_control$idx]

data_mean_control <- cbind(metadata_control,data.frame("intensities" = data_mean_control))

data_mean_control <- data_mean_control %>% dplyr::group_by(cell) %>% dplyr::summarize(mean_ion = median(intensities))

data_mean_drug <- data[iot,metadata_doi$idx]

data_mean_drug <- cbind(metadata_doi,data.frame("intensities" = data_mean_drug))

data_mean_drug <- data_mean_drug %>% dplyr::group_by(cell) %>% dplyr::summarize(mean_ion = median(intensities))

comb_data <- merge(data_mean_control, data_mean_drug, by = 'cell')

comb_data$perfect_fit <- comb_data$mean_ion.x

(ggplot(comb_data, aes(x=log10(mean_ion.x),y=log10(mean_ion.y), label = cell))+
    geom_point(col = 'red')+
    ggrepel::geom_text_repel(cex = 3)+
    theme_bw()+
    geom_abline(slope = 1, intercept = 0, col = 'red',linetype = "dashed")+
    geom_smooth(method='lm', formula= y~x) -> plt)



# R/S/I clustergram -------------------------------------------------------

# import R/S/I data
setwd(path_data_file)

RS_output <- read.csv("outcomes_GR24_RSgroups_filtered.csv")

#import FC data

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

lapply(unique(metab_fcs$drug), function(drug_idx){
  tmp <- subset(metab_fcs, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> metab_fcs

metab_fcs <- do.call(rbind, metab_fcs)

lapply(unique(RS_output$Drug), function(drug_idx){
  #drug_idx = "Erlotinib"
  RS_sub <- subset(RS_output, Drug == drug_idx & !is.na(percent_change_GR))

  for(groups_idx in list("R",c("S","I"))){

    #groups_idx <- c('S','I')
    groups_idx <- c('R')
    RS_sub_group <- subset(RS_sub, group %in% groups_idx)

    data_sub <- subset(metab_fcs, cell_line %in% RS_sub_group$cell & conc %in% RS_sub_group$Final_conc_uM)

    fcs_sub <- subset(metab_fcs, drug == drug_idx & conc %in% RS_sub_group$Final_conc_uM & cell_line %in% RS_sub_group$cell)

    #remove lowest concentration

    fcs_sub <- subset(fcs_sub, conc != 1)

    #prepare data as wide matrix format for clustering

    my_palette <- colorRampPalette(c("red", "white", "green"))(n = 299)

    fcs_sub$cell_conc <- paste(fcs_sub$cell_line, fcs_sub$conc, sep = "_")

    fcs_sub$log2fc <- ifelse(fcs_sub$log2fc < -2, -2, fcs_sub$log2fc)
    fcs_sub$log2fc <- ifelse(fcs_sub$log2fc >  2,  2, fcs_sub$log2fc)

    fcs_sub <- subset(fcs_sub, pvalue < 0.05)

    fcs_sub <- fcs_sub %>% group_by(ionIndex) %>%  filter(n() >= 6)

    fcs_sub <- tidyr::pivot_wider(fcs_sub, names_from = cell_conc, values_from = log2fc,id_cols = ionIndex,values_fill = 0)

    tmp_ions <- ions

    tmp_ions <- tmp_ions %>% dplyr::group_by(ionIndex) %>% slice(1)

    my_row_names <- tmp_ions[as.numeric(fcs_sub$ionIndex), "name"]$name

    fcs_sub$ionIndex <-NULL

    rownames(fcs_sub) <- my_row_names

    heatmaply::heatmaply(as.matrix(fcs_sub),
              density.info="none",  # turns off density plot inside color legend
              trace="none",          # turns off trace lines inside the heat map
              col = my_palette,na.rm = T,
              show_dendrogram = c(F, F),
              limits = c(-2,2))

    heatmaply::heatmaply(as.matrix(fcs_sub),
                         density.info="none",  # turns off density plot inside color legend
                         trace="none",          # turns off trace lines inside the heat map
                         col = my_palette,na.rm = T,
                         file = paste(drug_idx, groups_idx,'.html', sep = "_"))


    heatm(as.matrix(fcs_sub),
                      density.info="none",  # turns off density plot inside color legend
                      trace="none",          # turns off trace lines inside the heat map
                      col = my_palette,na.rm = T
                      )


  }







})


# export log2fc data for cytoscape ----------------------------------------

tmp <- metab_fcs

tmp <- subset(tmp, drug == 'Methotrexate' & cell_line == "SW620" & conc == 5, select = c('ionIndex', "log2fc",'pvalue'))

ions_sub <- ions[,c('ionIndex','idKEGG')]

ions_sub <- tidyr::separate_rows(ions_sub,idKEGG, convert = TRUE, sep = ' ')

ions_sub <- subset(ions_sub,!grepl("^\\s*$", idKEGG))

ions_sub <- ions_sub %>% dplyr::group_by(idKEGG) %>% dplyr::slice(1)

tmp <- merge(tmp, ions_sub,by='ionIndex')

write.csv(file = 'example_for_metabscape.csv', x = tmp[,c('idKEGG', 'log2fc','pvalue')])




# RS groups enrichment ---------------------------------------------------------------

# #calculate pathway enrichment for R/S groups ------------------------------------

lapply(unique(RS_groups$Drug), function(drug_idx){
  #drug_idx = 'BPTES'
  print(drug_idx)
  RS_sub <- subset(RS_groups, Drug == drug_idx & !is.na(percent_change_GR))

  res_sub <- subset(RS_sub, group == 'R' )

  sens_sub <-subset(RS_sub, group == 'S')

  int_sub <-subset(RS_sub, group == 'I')

  if(nrow(res_sub)==0 | nrow(sens_sub)==0){return(NULL)}

  if(nrow(res_sub) < nrow(sens_sub)){
    #add "I" to balance cases

    row_diff <- nrow(sens_sub)- nrow(res_sub)
    res_sub <- rbind(res_sub, int_sub[1:row_diff,])
    res_sub <- res_sub[!is.na(res_sub$X),]
    res_sub$group <- "R"
  }

  if(nrow(res_sub) > nrow(sens_sub)){
    #add "I" to balance cases

    row_diff <- nrow(res_sub)- nrow(sens_sub)
    sens_sub <- rbind(sens_sub, int_sub[1:row_diff,])
    sens_sub <- sens_sub[!is.na(sens_sub$X),]
    sens_sub$group <- "S"
  }

  meta_sub <- subset(metadata, drug == drug_idx)

  #TODO apply this code to the previous definition of conc_idxs
  meta_sub$conc <- meta_sub %>% dplyr::group_by(conc) %>%  group_indices(conc)
  #unique(meta_sub$tmp)
  #xtabs(~conc+tmp, meta_sub)


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

    ions_sub <- ions_sub %>% dplyr::group_by(idKEGG) %>% dplyr::slice(1)

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
    return(Output_GSEA)
  }


}) -> Output_GSEA


#plotting enrichment results

#setwd("C:\\Users\\masierom\\polybox\\Programing\\GSEA")
#files <- list.files(pattern = "SUMMARY.RESULTS.REPORT")

#results <- lapply(files, read.delim)

results <- Output_GSEA

names(results) <- unique(RS_groups$Drug)

lapply(seq_along(results), function(idx,x,y){

  x = unlist(x, recursive = F)
  named_result <- data.frame(do.call(rbind,x[idx]))
  if(nrow(named_result)>0){
    named_result$drug <- y[idx]
    return(named_result)
  }
},x = results, y = names(results)) -> results


results <- do.call(rbind, results)

results <- subset(results, glob.p.val < 0.05)

path_names <- lapply(results$GS,function(idx){
  KEGGREST::keggGet(idx)[[1]][[2]]})
results$path_names <- unlist(path_names)

results$path_names <- gsub(" - Homo sapiens (human)","",results$path_names,fixed = T)



#plotting results PEA

ggplot(results, aes(x= drug, y=path_names))+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))



# Compare metabolomics results for R/Sgroups, and see which one is better

# plot the metabolites that are increased in R or S as a volcano plot



# check which metabolites are driving the correlation, hopefully positive controls
# Correlate FC with GR50 and GR24
# For metabolic drugs: check if direct substrate/products involved in MoA are relating to any drug metric
# Multiomics: Link protein/gene/mRNA levels to metabolomics/drug sensitivity
# Multiomics: distance_to_target analysis
