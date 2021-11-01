## import data ####

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves.RData")

## import functions ####

source("C:\\Users\\masierom\\polybox\\Programing\\Project_exometabolites\\modelling_growth_curves.R")

library(dplyr); library(tidyr)

## get CV of growth rates at 24h post ttm ####

data_24 <- subset(data_corrected, Time %in% c(0,24))

data_24 <- subset(data_24, !Drug %in% c("exception"))

#data_24 <- subset(data_24, Drug == "DMSO")

group_names <- c("Drug", "Final_conc_uM", "cell")

groups <- data_24[,group_names]

groups <- groups %>% group_by_all() %>% slice(1) %>% ungroup()

rownames(groups) <- do.call(paste, groups)

groups$cv <- unlist(

lapply(rownames(groups), function(idx){
  
  #idx <- rownames(groups)[15]
  
  tmp_group <- groups[rownames(groups) == idx,]
  
  tmp_data <- subset(data_24, 
                     Drug == tmp_group$Drug &
                     Final_conc_uM == tmp_group$Final_conc_uM &
                     cell == tmp_group$cell)
  
  diff_conf <- lapply(unique(tmp_data$Well), function(well){
    
    #well = unique(tmp_data$Well)[3]
    
    tmp_data_well <- tmp_data[tmp_data$Well == well,]
    
    diff <- subset(tmp_data_well, Time == max(tmp_data_well$Time), Conf, drop = T) - 
      subset(tmp_data_well, Time == min(tmp_data_well$Time), Conf, drop = T)
    
    return(diff)
  })
  
  diff_conf <-  unlist(diff_conf)
  
  sd_group <- sd(diff_conf)/mean(diff_conf)
  
  return(sd_group)
  
}))

groups$cv <- abs(groups$cv) * 100

groups <- groups[groups$cv < 100,]

hist(groups$cv)

sum(groups$cv <20, na.rm = T) / nrow(groups) * 100

library(ggplot2)

ggplot(subset(groups, Drug == "Docetaxel"), aes(cv, col = factor(Final_conc_uM)))+
  geom_density()+
  theme_bw()+
  theme(legend.position = "none")
  # facet_wrap(~cell)

ggplot(subset(groups), aes(x = log(cv), y = log(Final_conc_uM)))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  facet_grid(~Drug)


## get growth rate +/- 30 minutes around sampling point ####

lapply()


### plot stuff #####

cc <- scales::seq_gradient_pal("lightblue", "red")(seq(0,1,length.out= 6))


tmp_fig <- subset(data_corrected, Time <= 72 & Drug %in% c("Oligomycin A", "DMSO"))

tmp_fig$Final_conc_uM <- ifelse(tmp_fig$Final_conc_uM == 367, "DMSO", tmp_fig$Final_conc_uM)

tmp_fig$Final_conc_uM <- factor(tmp_fig$Final_conc_uM)

names(tmp_fig)[3] <- "Confluence"

ggplot(tmp_fig,
       aes(x = Time, y = Confluence, color = factor(Final_conc_uM), linetype = Drug, group = factor(Final_conc_uM)))+
  geom_smooth(size = 1.3)+
  #scale_color_manual(values = cc)+
  #facet_wrap(~Drug+cell)+
  theme_bw()+
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())+
  theme(panel.border     = element_blank())+
  theme(axis.ticks       = element_blank())+
  theme(panel.border     = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=20))+
  facet_wrap(~cell)+
  geom_vline(xintercept = 48, col = "red")+
  geom_vline(xintercept = 24, col = "blue")

ggsave(p1)

rm(tmp_fig)

data_corrected$Time <- round(data_corrected$Time,0)

ggplot(subset(data_corrected, source_plate == "P1.txt" & Time <= 72 & Drug %in% c("DMSO")),
       aes(x = Time, y = Conf, color = factor(Well), group = factor(Well)))+
  geom_line()+
  facet_grid(~cell) #FIXME When we plot DMSO by cell line, we have a few outliers. Remove these outliers.

ggplot(subset(data_corrected, Time <= 72 & Drug %in% c("Pemetrexed", "DMSO") & cell == "M14"),
       aes(x = Time, y = Conf, color = factor(Well), group = factor(Well)))+
  geom_line()+
  facet_wrap(~ Final_conc_uM) #FIXME When we plot one cell line, we have a few outliers. Remove these outliers.

#TODO plot variation across replicates by well in the 384 well plate, hopefully only certain wells in the edge will have problems

# maybe calculate intra replicate CV adn use ANOVA and check the variance by well, see if its consistant across CLs or if it
#varies from CL to CL (plate to plate). Based on this, either exclude wells or 

#TODO based on the plots above, create an outlier removal function that will remove the obvious outliers

# https://www.r-bloggers.com/outlier-detection-and-treatment-with-r/
# https://stepupanalytics.com/outlier-detection-techniques-using-r/ 

# Calculate Metrics of Growth Inhibition --------------------------------
# Calculate GI50, IC50, GRmetrics 50
