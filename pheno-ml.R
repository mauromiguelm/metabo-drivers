setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data")

library(ggplot2)
library(dplyr)

source("C:\\Users\\masierom\\Downloads\\tdsr-pheno-ml\\core.R")
source("C:\\Users\\masierom\\Downloads\\tdsr-pheno-ml\\fileConversion.R")

phenoDist <- "ACHN_CL3_P1.json"

data <- phenoMLtoCaseB(phenoDist)

data$concentration <- round(as.numeric(data$concentration), 5)

data$cell_line <- strsplit(phenoDist, split = "_")[[1]][1]

output <- tdsR::tdsR_fit(inputData = data, groupingVariables = c("agent", "cell_line"), case = "B", timeTreatment = 10, smoothData = T,
               lowerLimit = 0.05, limitThreshold = 0.2, orderConc = T)

output_tds <- tdsR::tdsR_getOutput(output, "tdsR")

data$groups <- do.call(paste, data[,c("agent", "cell_line")])

data_comb <- data %>% dplyr::left_join(output_tds, by = "groups")

ggplot(data_comb, aes(x = time, y= value, col = factor(concentration)))+
  geom_smooth()+
  geom_vline(aes(xintercept = tds))+
  facet_wrap(~agent, ncol = 3)+
  scale_colour_hue(h = c(270, 360))+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  geom_vline(xintercept = 0, linetype="dotted", size = 1)
