# load libraries, paths ----------------------------------------------------------

path_data_file = "C:\\Users\\mauro\\Documents\\phd_results"
path_fig = "C:\\Users\\mauro\\Documents\\phd_results"

# import clustering results -----------------------------------------------

files <- list.files(paste0(path_data_file, "\\cluster_labels"), full.names = T)

cluster_data <- lapply(files, read.csv)

# import omics ------------------------------------------------------------

setwd(paste(path_data_file, 'omics',sep = "//"))

omics_type <- list.files()

omics_data <- lapply(omics_type,read.csv)

discrete_omics <- c(1,3,5)
  
continuous_omics <- c(2,4,6,7)


# run associations --------------------------------------------------------

#discrete predictors

lapply(1:length(cluster_data), function(drug_idx){
  
  results <- data.frame()
  #print(paste("drug =",drug, "drug_index =", drug_idx))
  #drug_idx = 1
  #TODO iterate over concentration
  
  cluster_of_drug <-  cluster_data[[drug_idx]]
  drug <- unique(cluster_of_drug$drug)
  cluster_of_drug <- cluster_of_drug[cluster_of_drug$conc==5,] #only check clusters at high concentration
  cluster_of_drug <- cluster_of_drug[,c("cell_line", "clusters")]
  cluster_of_drug <- cluster_of_drug[cluster_of_drug$clusters !=-1,]
  
  if(length(unique(cluster_of_drug$clusters))>1){
    #only proceed if max concentration has >1 cluster label
    for(omic in continuous_omics){
      #omic = 2
      #print(paste("omic =", omic))
      data <- omics_data[[omic]]
      colnames(data)[2] <- "cell_line"
      data <- merge(cluster_of_drug, data, by ="cell_line")
      feature_names <- colnames(data)
        
      for(feature in 8:dim(data)[2]){
        #print(paste("feature =", feature))
        #feature = 8
        data_feature <- data[,c(1:7, feature)]
        
        if(sum(!is.na(data_feature[,8]))>6){
          #print(paste("feature =", feature, "_passed"))
          #only proceed if binary feature has >1 label
          
          stat_test <-  try(aov(data_feature[,8]~data_feature$clusters), silent = T)
          p_val <- summary(stat_test)[[1]]
          p_val <- p_val$`Pr(>F)`[1]
          results <- rbind(results, data.frame(drug,omic, feature,name=feature_names[feature],  p_val))
          
          } 
      }
    }
  }
  
return(results)
}) ->results_discrete

tmp <- do.call(rbind, results_discrete)

tmp$p_val <- ifelse(tmp$p_val==0,min(tmp$p_val[tmp$p_val!=0], na.rm = T)*0.1,tmp$p_val )

tmp$p_val_adj <- p.adjust(tmp$p_val, method = 'fdr')

tmp$pass <- ifelse(tmp$p_val_adj<0.05, 1,0)

tmp$omic_feature <- paste(tmp$omic, tmp$feature)

aggregate(pass~omic,tmp, sum, na.rm = T)

aggregate(pass~drug,tmp, sum, na.rm = T)

tmp1 <- aggregate(pass~omic_feature,tmp, sum, na.rm = T)

sum(tmp$pass, na.rm = T)

write.csv(subset(tmp, pass ==1), 'associations_omics.csv')

write.csv(tmp, 'associations_omics_full.csv')

# define ranking

#tmp <- tmp[!is.na(tmp$p_val),]
#tmp <- tmp[!tmp$p_val==1,]

library(tidyr)
library(dplyr)
lapply(unique(tmp$omic), function(x){
  #x=2
  tmp_omic <- subset(tmp, omic ==x)
  tmp_omic <- tmp_omic %>% mutate(rank = dense_rank(p_val))
  return(tmp_omic)
  
}) -> tmp

tmp <- do.call(rbind, tmp)

# plot associations by drug

library(ggplot2)



#ggplot(tmp, aes(rank,p_val, col = omic))+
#  geom_point()

pdf("omics_ranking.pdf",width = 4, height = 2)
ggplot(subset(tmp, omic %in% c(2,6)), aes(log10(rank),log10(p_val), col = factor(omic)))+
  geom_line()+
  theme_bw()

dev.off()


