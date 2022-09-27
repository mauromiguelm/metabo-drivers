## This code analyse metabolomics data

# load packages and definitions -------------------------------------------

library(openxlsx); library(dplyr); library(rhdf5); library(ggplot2); library(parallel); 
library(proxy); library(KEGGREST); library(viridis); library(circlize); library(confluencer)
library(ggrepel);

#import cell metadata from depmap 

cell_metadata <- read.csv("metadata/metadata_cells_depmap.csv")

colnames(cell_metadata)[2] = 'cell'

#import drug sensitity metrics

data_GR50 <- read.csv("data/confluence/outcomes_growth_inhibition50.csv")

#import kegg hsa compounds

kegg_hsa_cpds = read.csv("metadata/kegg_hsa_cpds.csv")

#import metabolomics

data<- rhdf5::h5read(file = "data/metabolomics/metabolomics.h5", '/data')

metadata <- rhdf5::h5read(file = "data/metabolomics/metabolomics.h5", '/metadata')

metadata$X <- NULL

metadata$idx <- as.numeric(metadata$idx)

metadata <- metadata[order(metadata$idx),]

rownames(metadata) <- NULL

ions <- rhdf5::h5read(file = "data/metabolomics/metabolomics.h5", '/ions')

ions <- data.frame(ions)

# impot R/S groups

RS_groups <- read.csv("data/outcomes_GR24_RSgroups_filtered.csv"  )

#define drugs in screen

drugs_in_screen <- c(unique(metadata$drug)[!(unique(metadata$drug) %in% c("PBS", "DMSO","hct15Ctrl", "SolvCrtl","hct15Mtx","poolCtrl",'NA'))])

# plot a metabolomics expected patterns in controls -------------------------------------

# glucose vs. lactate

hexoseP_idx = which(ions$name=='Hexose')

lactate_idx = which(ions$name=='Lactate')

dmso_controls = metadata[which(metadata$drug=='DMSO' & metadata$conc=='367'),'idx']

data_control = data.frame(t(data[c(hexoseP_idx, lactate_idx),dmso_controls]))

data_control$idx = dmso_controls

colnames(data_control)[1:2] = c("hexose", 'lactate')

data_control = merge(data_control, metadata, by="idx")

data_control = data_control[,c('cell',"hexose", 'lactate')] %>% group_by(cell) %>% summarise_all(list(mean, sd))

data_control <- merge(data_control,cell_metadata, by = 'cell')

fit <- lm(scale(lactate_fn1)~scale(hexose_fn1),data_control)

data_control$colors_lineage1 <- factor(data_control$lineage_1,labels = colorspace::diverge_hcl(length(unique(data_control$lineage_1))))

coefs <- coef(fit)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x))

pdf('figures/hexose_lactate_fit.pdf',width = 3,height = 3.5)

plot(x = scale(data_control$hexose_fn1), y= scale(data_control$lactate_fn1),
     xlab='',ylab = '',col = data_control$colors_lineage1, pch=19,ylim = c(-2,2.5),xlim = c(-2.5,2))

col_df = data_control%>% group_by(lineage_1) %>% slice(1)

abline(fit,col='black',lty=1,lwd=1)

legend('topleft',title="Lineage", cex=0.1, pch=19,
       col=col_df$colors_lineage1,
       legend=col_df$lineage_1, ncol=1)

text(-2.5, 2.0, eqn, pos = 4)

dev.off()

# calculating log2(FCs) for all data and for every drug.  --------

#iterating over plates, cells, drug, and conc to calculate FCs

metadata$cell_plate <- paste(metadata$cell, metadata$source_plate, sep = "_")

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

    write.csv(metab_fcs, paste0("data//metabolomics//log2fc//",cell_plate_idx,'.csv'))

  }else{
    NA
  }
})

# export log2fc data for clustering -------------------------------------

#map concentrations

metab_fcs <- lapply(list.files("data//metabolomics//log2fc",pattern = "_P"),read.csv)

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

groups_to_keep <-read.csv('data//outcomes_GR24.csv')

groups_to_keep$cell_drug_conc <- paste(groups_to_keep$cell,groups_to_keep$Drug,groups_to_keep$Final_conc_uM)

#metab_fcs <- subset(metab_fcs, cell_drug_conc %in% unique(groups_to_keep$cell_drug_conc))

#return a matrix of data (rows are samples, cols are ion features) and metadata (sample, growth inhibition)

lapply(unique(metab_fcs$drug), function(drug_idx){
  #drug_idx = "Methotrexate"
  print(drug_idx)
  tmp <- subset(metab_fcs, drug == drug_idx)

  tmp$pvalue <- NULL

  tmp <- tidyr::spread(tmp, key = ionIndex,value = log2fc)

  metadata <- tmp[,1:6]

  data <- tmp[,7:ncol(tmp)]

  metadata <- merge(metadata, groups_to_keep[,c('cell_drug_conc','percent_change_GR')], by ='cell_drug_conc',all.x = T)

  write.csv(metadata, paste("data//metabolomics/log2fc//metadata",drug_idx,"log2fc.csv",sep ="_"))
  write.csv(data, paste("data//metabolomics/log2fc//data",drug_idx,"log2fc.csv",sep ="_"))
})


# association GR24 with log2fc --------
# lm(metab~GR24) across all cell lines and concentrations that we have filtered strong effects/unnefective concentrations
# import log2fc data

metab_fcs <- lapply(list.files(path = 'data//metabolomics//log2fc',pattern = "_P"),read.csv)

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
  #drug_idx = 'Methotrexate'

  #distribute nboot across all ions
  print(paste('drug:',paste(drug_idx),'is now running'))
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

    perm_df <- t(replicate(max(iter_per_ion), sample(length(unique(tmp$cell_conc)), replace = FALSE, prob = NULL)))

    #for every ion in each drug associate with GR50

    for(ion_idx in 1:nions){
      #ion_idx=3
      print(paste(drug_idx,ion_idx))

      tmp_ion <- subset(tmp, ionIndex == ion_idx)

      count_iter <- 0

      while(count_iter < iter_per_ion[ion_idx]){
        print(count_iter)

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

names(results) <- unique(metab_fcs$drug)

save(results,file = 'data//metabolomics//log2fc//bootstrap_results.Rdata')

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

write.csv(bootstrap_ci,file = 'data//metabolomics//log2fc//confidence_intervals_from_permutation_1.csv')

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

write.csv(slope_metabolite_effect_on_growth, 'data//metabolomics//log2fc//metabolite_GR24_association.csv')

#select the ions that survived permutation thresohld

df <- read.csv('data//metabolomics//log2fc//metabolite_GR24_association.csv')

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

write.csv(df, 'data//metabolite_GR24_association_survivingCI.csv')

# baseline vs drug treated ion compairison --------------------------------

lapply(1:nrow(ions), function(ionidx){
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

write.csv(basal_associations, 'data//basal_association.csv')

# Determining half effect concentration for relevant ions -----------------

#import GR24 associations with metabolism

#read associations

df <- read.csv('data//metabolite_GR24_association_survivingCI.csv')

df$X <- NULL

#import full GR24 data

outcomes_GR24 <- read.csv("data//outcomes_GR24.csv")

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

write.csv(drug_cell_fc, "data//drug_cell_fc.csv")

write.csv(drug_cell_ED50,"data//drug_cell_ED50.csv")

#TODO save models, so I can use them later

#read results

drug_cell_fc <- read.csv("data//drug_cell_fc.csv")

drug_cell_ED50 <- read.csv("data//drug_cell_ED50.csv")

#gap fill NA with high concentration

lapply(unique(drug_cell_fc$idx), function(idx){
  #idx = 10
  sub_data <- drug_cell_fc[drug_cell_fc$idx==idx,]
  
  sub_data$ED50 <- ifelse(is.na(sub_data$ED50), max(sub_data$ED50, na.rm = T)*1.1, sub_data$ED50)
  
  return(sub_data)
}) -> tmp

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
    ggsave(filename = paste0("figures//",sub_fc$drug,metab_name,".png"), plot = plot)
    
  }
  
  return(data.frame(unique(sub_data$idx),GR_diff))
  
}) -> tmp

EC50_statistics <- data.frame(do.call(rbind,tmp))

colnames(EC50_statistics) <- c("association_index","EC50_FC")

write.csv(EC50_statistics, "data\\statistics_EC50.csv")

# combine all metrics into one matrix -------------------------------------

#read EC50 results

drug_cell_fc <- read.csv("data\\drug_cell_fc.csv")

drug_cell_ED50 <- read.csv("data\\drug_cell_ED50.csv")

EC50_statistics <- read.csv("data\\statistics_EC50.csv")

EC50_statistics$X <- NULL

#select the ions that survived permutation thresohld

df <- read.csv('data\\metabolite_GR24_association_survivingCI.csv')

df$X <- NULL

df$color <- ifelse(df$slope<0, "Negative", "Positive")

colnames(df)[4:8] <- paste("DMA", colnames(df)[4:8], sep = "_") 

df <- merge(df, ions, by='ionIndex')

write.csv(df, file = 'data\\metabolite_GR24_association_survivingCI_ionLabel.csv')

# combine log2fc/GI50 association with basal/treated metabolome association

df$association_index = df$idx

df$drugion <- paste(df$drug, df$ionIndex)

basal_associations <- read.csv('data\\basal_association.csv')

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

write.csv(df, "data\\drug_metabolite_associations.csv")

#generate statistics as results per group

df <- read.csv("data\\drug_metabolite_associations.csv")

stats_associations <-  xtabs(~Basal_effect+DMA+GR_independence, df)

#generate number of associations per ion

stats_count_per_ion <- df %>% group_by(name) %>% summarize(count = n())


# plot associations by drug, heatmap by metrics ---------------------------

#import results

df <- read.csv("data\\drug_metabolite_associations.csv")
df$X = NULL

lapply(unique(df$drug), function(drug){
  #drug = "Decitabine"
  tmp <- df[df$drug==drug,]
  
  n_interactions <- length(unique(tmp$name))
  
  #subset ion and the tree metrics
  
  
  tmp <- tmp[,c("name","mzLabel","EC50_FC","DMA_slope","basal_slope")]
  
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
  
  ggsave(paste0("figures\\",drug, "_h1.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.34,width = 7, limitsize = F)
  ggsave(paste0("figures\\",drug, "_h0.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.4,width = 7, limitsize = F)
  ggsave(paste0("figures\\",drug, "_h2.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.5,width = 7, limitsize = F)
  ggsave(paste0("figures\\",drug, "_h3.pdf"),plot = plt,height = length(unique(tmp$Metabolite))*0.3,width = 9, limitsize = F)
  
  #create a legend for GR independence
  
  GRI <- seq(c(-GRI_max), GRI_max, length.out = 7)
  
  GRI_col <- GRI_colmap(GRI)

  GRI <- round(GRI,0)
  
  pdf(paste0("data\\",drug, "_legend.pdf")) 
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x = 0, y = 1, fill = GRI_col, legend = GRI, title = "log2(GR independence)")
  dev.off()

  })


# plot top ranked interesting associations -----------------------

#number of associations per metabolite with diferent drugs (plot)

tmp <- df %>% group_by(ionIndex) %>% dplyr::summarise(count_metab = n())

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

ggsave("figures//number_metabs_associated_with_number_drugs.pdf",plt,device = "pdf",width=5, height=4,scale = 0.9)

#clirclize plot of top 10 associations ---------------------------------

df <- read.csv("data\\drug_metabolite_associations.csv")

df <- df[ifelse(df$name =="",F,T),]  #remove non annotated

drugs_to_plot <- df %>% group_by(drug) %>% filter(n()>=5)

df <- df[df$drug %in% drugs_to_plot$drug,] #remove drugs with less than 5 associations

#get the top 10 DMA associations

df <- df %>% group_by(drug) %>% dplyr::arrange(abs(DMA_slope)) %>% dplyr::slice_tail(n=10)

df <- df %>% group_by(drug) %>%dplyr::arrange(abs(DMA_slope)) %>%  dplyr::mutate(slope_norm = c(10:1)[1:n()])

#define boundaries

df$basal_slope <- ifelse(df$basal_slope >=2,2, df$basal_slope)
df$basal_slope <- ifelse(df$basal_slope <=c(0),-0, df$basal_slope)

df$DMA_slope <- ifelse(df$DMA_slope >=2,2, df$DMA_slope)
df$DMA_slope <- ifelse(df$DMA_slope <=c(-2),-2, df$DMA_slope)

#shorten the ion name so it fits in plotting

df[which(df$mzLabel == "mz579.0267"),'name'] <- "UDP glucuronic acid" #Uridine diphosphate glucuronic acid

df[which(df$mzLabel == "mz377.0325"),'name'] <- "Me-TIMP" #6-Methylthiopurine 5''-monophosphate ribonucleotide

df[which(df$mzLabel == "mz363.0170"),'name'] <- "6-Thio-IMP" #6-Thioinosine-5''-monophosphate

df[which(df$mzLabel == "mz613.1403"),'name'] <- "CMP-NeuNAc" #Cytidine monophosphate N-acetylneuraminic acid

df$name <- paste(df$name, paste0("(",df$mzLabel,")"), sep = " ")

df$EC50_FC <- log2(df$EC50_FC)

df$EC50_FC <- ifelse(df$EC50_FC >4,4, df$EC50_FC)
df$EC50_FC <- ifelse(df$EC50_FC <c(-4),-4, df$EC50_FC)

# merge MoA 

drug_metadata <- read.csv("data\\drug_metadata_moa.csv")

colnames(drug_metadata)[5] <- 'drug'

df <- merge(df, drug_metadata[,c('drug', 'Class')], by = 'drug')

#start circos plot

df$x=0

df$color <- ifelse(df$DMA_slope<0, "#FF0000", "green4")

png("figures//association_GR24_metabolites.png",width = 7000,height = 7000,res = 700)

x1 <- x2 <- y1 <-y2 <-0.5

circos.par("track.height" = 0.05,canvas.xlim = c(-(1+x1), 1+x2), canvas.ylim = c(-(1+y1), 1+y2),
           start.degree = 93)
circos.initialize(df$drug, x = (df$slope_norm))

circos.track(df$drug, y =df$DMA_slope,
             panel.fun = function(x, y) {
               #circos.axis(major.tick = F,labels = NULL)
             },bg.border = NA)

circos.trackText(sectors = df$drug,x = df$slope_norm, y= df$x- mm_y(3), labels = df$name,
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
                           CELL_META$cell.ylim[2] - mm_y(2),
                           CELL_META$sector.index,cex = 1,facing = 'clockwise',niceFacing = T, adj = c(1, 0.5))
             },bg.border = NA)

dev.off()
circos.clear()



col3 = colorRamp2(c(-GRI_max, 0, GRI_max), c("darkgreen", "grey", viridis(100)[100]))

df$color3 <- col3(df$EC50_FC)


pdf(paste("figures\\circle-legend-GRI.pdf")) 
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = 0, y = 1, fill = col3(seq(-GRI_max,GRI_max,1)),legend = seq(-GRI_max,GRI_max,1), title = "log2(GR independence)")
dev.off()


pdf(paste0("figures\\circle-legend-basal.pdf")) 
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

ggsave(paste0("figures\\", "fraction_GI50.pdf"),plot = plt,  width = 25, height = 15, units = "cm")
  
# plot summary all associations in one plot -----------------------------

df$drugmetab <- paste(df$drug, df$name)

df$drugmetab = ifelse(df$drug %in% c('Decitabine', 'BPTES'), df$drugmetab,"")

df$drugmetab = ifelse(abs(df$DMA_slope)>2.5, df$drugmetab,"")


ggplot(df, aes(x = rank(DMA_slope) ,y=DMA_slope, label = drugmetab))+
  geom_point()+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf,size = 3)+
  theme_bw()+
  xlab("Ranked DMA")+
  ylab("DMA")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))->plt


ggsave(filename = 'figures//rankedDMA.pdf',plt,width =3,height = 2.6)
