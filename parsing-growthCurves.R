
library(ggplot2); library(GRmetrics); library(RColorBrewer)

# import functions

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\CellCultureAnalysis.R")

source("C:\\Users\\masierom\\polybox\\Programing\\Project_exometabolites\\modelling_growth_curves.R")

# Importing source plates

source("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cell_culture_data\\190310_pre_screen_test_cpddilution\\exceptions\\log_processing.r")

source_plates <- data.frame(filenames = list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cell_culture_data\\190310_pre_screen_test_cpddilution\\20190513\\results", pattern = ".txt"))

source_plates$sourceid <- ifelse(grepl(pattern = "P1", source_plates$filenames), "1MSP001", "2MSP001")

#import growth curves

setwd("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_Pre_screen_test_CPDdilution\\20190513\\results")


fileNames <- list.files(full.names = F, pattern = ".txt")

data <- 

lapply(fileNames, function(x) {
  
  df <-
  ReadIncuCyteData(read_platemap = F, FileName_IncuCyte = x, Plate_size = 384, FileDirectory = getwd())
  
  df$Time <- as.numeric(df$Time)
  
  return(df)
  
  })

names(data) <- fileNames


# correcting time for HCT-15: HCT-15 did not record some wells, so I am correcting for this

data_corrected <- data

data_corrected$HCT15_P1.txt$Time <- data_corrected$HCT15_P1.txt$Time +
  (max(data_corrected$A549_P1.txt$Time) - max(data_corrected$HCT15_P1.txt$Time))


# removing bad wells based on first 24h  growth

store_HCT15 <- data_corrected$HCT15_P1.txt # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data

data_corrected <-
  
  lapply(names(data_corrected), FUN = function(x){
    
    #x = "A549_P1.txt"
    
    plate_name = x
    
    x = data_corrected[[x]]
    
    good_wells <-
      filter_growth_outliers(plate_name = plate_name, data = x, time_control = 24, save_diag_plots = T,
                             save_plots_directory = "\\\\imsbnas.d.ethz.ch/sauer1/users/mauro/cell_culture_data/190310_pre_screen_test_cpddilution/figures")
    
    new_data <-
      x[!(x[,"Well"] %in% good_wells$wells_exception),]
    
    print(paste("Plate:", plate_name,
                "Exception-Number:",good_wells$wells_exception_number
                ))
    
    return(new_data)
    
  })

names(data_corrected) <- fileNames

data_corrected$HCT15_P1.txt <- store_HCT15  # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data

rm(store_HCT15)

# Fit the confluence for each well, and return fitted confluence. Keep the same plate map structure.

base::lapply(names(data_corrected), function(plate_name){
  
  #plate_name =  names(data_corrected)[3] #FIXME delete me
  
  r.2.threshold = 0.8
  
  data_raw <- data_corrected[[plate_name]]
  
  fitted_data <- 
  
  lapply(unique(data_raw$Well), function(idx_well){
    
    #idx_well <- unique(data_raw$Well)[2] #FIXME delete me 
    
    data_well_raw <- subset(data_raw[,c("Time", "Conf", "Well")], Well == idx_well)
    
    min_scan_Time <- ceiling(min(data_well_raw$Time))
    
    max_scan_Time <- floor(max(data_well_raw$Time))
    
    model_pred_poly <- get_growthMetrics(data_well_raw)[[3]]
    
    time_sequence <- seq(min_scan_Time ,max_scan_Time,1)
    
    model_well_pred <- predict(model_pred_poly, data.frame(Time = time_sequence))
    
    pred_conf <- data.frame(Well = idx_well, Time = time_sequence, Conf = model_well_pred)
    
    well_metrics <- summary(model_pred_poly)
    
    well_metrics <- data.frame(Well = idx_well, adj.r.2 = well_metrics$adj.r.squared)
    
    return(list(pred_conf, well_metrics))
    
  })
  
  # save diagnostic plots
  
  diagnostics_well <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[2]]) ) #for each plate, plot diagnostic plots with the R-squared.. check if any of fits have ploblems
  
  save_plots_directory = "\\\\imsbnas.d.ethz.ch/sauer1/users/mauro/cell_culture_data/190310_pre_screen_test_cpddilution/figures"
  
  filename = paste("fitting-QC", plate_name, gsub(" ", "_",as.character(Sys.time())), ".png", sep = "_")
  
  filename = gsub(pattern = ":", replacement = "_", x = filename)
  
  bad_wells = diagnostics_well$Well[diagnostics_well$adj.r.2 < r.2.threshold]
  
  #print(paste(as.character(bad_wells), plate_name), sep = "_")
  
  png(filename = paste(save_plots_directory, filename , sep = "/"))
  fig <- hist(diagnostics_well$adj.r.2, main = plate_name)
  text( paste0("Adj.r.2 < 0.8 = ", length(bad_wells)) , x = fig$breaks[3], y = max(fig$counts))
  dev.off()
  
  # return a df with the right columns
  
  pred_conf <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[1]]) ) #for each plate, plot diagnostic plots with the R-squared.. check if any of fits have ploblems
  
  pred_conf <- subset(pred_conf, !(Well %in% bad_wells))
  
  return(pred_conf)
  
}) -> data_corrected


names(data_corrected) <- fileNames

# combining metadata of source.plates into growth_data


source_layout <- # import source plate layout data
  
  lapply(list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cpd_data\\190508_msp_layout",
                    pattern = "randomized",
                    full.names = T),
         function(x) {
           
           #x = "\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cpd_data\\190508_msp_layout/randomized_layout_2MSP.xls"
           
           tmp_data  = readxl::read_xls(x)
           
           tmp_data$Well = paste0(tmp_data$Row, tmp_data$Column)
           
           return(tmp_data)
           
         })


data_corrected <-
  
  lapply(names(data_corrected), function(x){
    
    stopifnot(require(dplyr))
    
    tmp_data = data_corrected[[x]]
    
    #x = "A549_P1.txt"
    
    match_source <- source_plates[source_plates$filenames == x,"sourceid"] 
    
    if(grepl(pattern = "1MSP", x = match_source)){
      
      tmp_data <- inner_join( tmp_data, source_layout[[1]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
      
      return(tmp_data)
      
    }else if(grepl(pattern = "2MSP", x =match_source)){
      
      tmp_data <- inner_join(  tmp_data, source_layout[[2]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
      
      return(tmp_data)
      
    }else{
      
      stop("some experimental plate could not be matched to a source plate layout.")
      
    }
    
  } )

names(data_corrected) <- fileNames

# relabelling bad wells based on Echo transfer. I can relabel those with medium, to increase the number of controls

data_corrected <-
  
  lapply(names(data_corrected), function(x){
    
    #x = "A549_P1.txt"
    
    source_id <- unlist(subset(source_plates, filenames == x, sourceid))
    
    exception_wells <- subset(exceptions, Destination.Plate.Name  == source_id, Destination.Well, drop = T)
    
    tmp_data <- data_corrected[[x]]
    
    tmp_data[(tmp_data$Well %in% exception_wells), c("Drug")] <- c("exception")
    
    tmp_data[(tmp_data$Well %in% exception_wells), c("Final_conc_uM")] <- c(333)
    
    tmp_data$cell<- strsplit(x, split = "_")[[1]][1]
    
    tmp_data$source_plate <- strsplit(x, split = "_")[[1]][2]
    
    return(tmp_data)
    
  })

names(data_corrected) <- fileNames

# checking if exceptions show the same behavious as PBS or DMSO

data_corrected <- do.call(rbind, data_corrected)

# removing DMSO that is not 33.3 and 3.33, as it causes many data analysis problems

data_corrected <- subset(data_corrected, !(Drug == "DMSO" & Final_conc_uM != 333) )

# plotting all data to check for inconsistencies

#$cc <- scales::seq_gradient_pal("lightblue", "red")(seq(0,1,length.out= length(unique(data_corrected$Final_conc_uM))))

cc <- scales::seq_gradient_pal("lightblue", "red")(seq(0,1,length.out= 9))

ggplot(subset(data_corrected, source_plate == "P1.txt" & Time <= 72 & cell == "A549" & Drug %in% c("Docetaxel", "DMSO")),
       aes(x = Time, y = Conf, color = factor(Final_conc_uM), group = factor(Final_conc_uM)))+
  geom_smooth()+
  scale_color_manual(values = cc)+
  facet_wrap(~Drug+cell)

data_corrected$Time <- round(data_corrected$Time,0)  

ggplot(subset(data_corrected, source_plate == "P1.txt" & Time <= 72 & Drug %in% c("DMSO")),
       aes(x = Time, y = Conf, color = factor(Well), group = factor(Well)))+
  geom_line()+
  facet_grid(~cell) #FIXME When we plot DMSO by cell line, we have a few outliers. Remove these outliers.

ggplot(subset(data_corrected, source_plate == "P1.txt" & Time <= 72 & Drug %in% c("Pemetrexed", "DMSO") & cell == "NCIH460"),
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



# GRmetrics ---------------------------------------------------------------


# GR metrics, their example code only use the endpoint and time at treatment, but I will try all time points

# we need to select time point zero as the time we treated the samples (as per their requirements)


#time_treatment <- 24


drugs_in_screen <- "17-AAG"  #c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO", NA, "Water", "WATER", "H2O", "Dmso", "Control", "Ctrl", "exception"))]))

matrix()   # columns will be drugs, rows will be cell lines

data_GRmetrics <- data_corrected

cell_line <- "NCIH460"
drug  <- "Docetaxel"
time_treatment <- 24

data_Grmetrics_Ttm <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # preparing the data for one drug
data_Grmetrics_Ctr <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # getting the matching controls ready

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO" & Final_conc_uM == 333) # I separated the drug from the control, as its easier to parse

# keep only matching time points

match_time_intervals <- intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time) #make sure both datasets cover the same time

time_treatment <- max(match_time_intervals[match_time_intervals < time_treatment])

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)

# Keeping certain constraints on time, so we are sure conf is not limiting (eg. >80)
# also, GRmetrics require the time point at treatment, and not before. So zero will be time point slightly before treatment

data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time <= 72)
data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time <= 72)

data_Grmetrics_Ctr$Time <- data_Grmetrics_Ctr %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
data_Grmetrics_Ttm$Time <- data_Grmetrics_Ttm %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time



#data_Grmetrics_Ctr$Time <- ifelse(data_Grmetrics_Ctr$Time == min(data_Grmetrics_Ctr$Time), 0, data_Grmetrics_Ctr$Time) #FIXME for each cell line get a specific time, instead of general
#data_Grmetrics_Ttm$Time <- ifelse(data_Grmetrics_Ttm$Time == min(data_Grmetrics_Ttm$Time), 0, data_Grmetrics_Ttm$Time)

#TODO select only one time point, eg. 60h, and then calculate the GRmetric across this time point.. could be even 48h

# Prepare data_Grmetrics_Ctr as example Case C

data_Grmetrics_Ctr$time = data_Grmetrics_Ctr$Time

data_Grmetrics_Ctr$cell_line <- (do.call(rbind, strsplit(x = rownames(data_Grmetrics_Ctr), "_")))[,1]

data_Grmetrics_Ctr$agent = "-"

data_Grmetrics_Ctr$concentration = 0

data_Grmetrics_Ctr$cell_count <- data_Grmetrics_Ctr$Conf

data_Grmetrics_Ctr <- data_Grmetrics_Ctr[,c("cell_line", "agent","time", "concentration", "cell_count")]

rownames(data_Grmetrics_Ctr) <- NULL

# Prepare data_Grmetrics_Ttm as example Case C

head(data_Grmetrics_Ttm)

data_Grmetrics_Ttm$time = data_Grmetrics_Ttm$Time

data_Grmetrics_Ttm$cell_line <- (do.call(rbind, strsplit(x = rownames(data_Grmetrics_Ttm), "_")))[,1]

data_Grmetrics_Ttm$agent <- (data_Grmetrics_Ttm$Drug)

data_Grmetrics_Ttm$concentration = data_Grmetrics_Ttm$Final_conc_uM

data_Grmetrics_Ttm$cell_count <- data_Grmetrics_Ttm$Conf

data_Grmetrics_Ttm <- data_Grmetrics_Ttm[,c("cell_line", "agent","time", "concentration", "cell_count")]

rownames(data_Grmetrics_Ttm) <- NULL

# merge drug df and control df

data_comb <- rbind(data_Grmetrics_Ctr, data_Grmetrics_Ttm)

output1 = GRfit(inputData = data_comb, groupingVariables = 
                  c("cell_line", 'agent', 'time'), case = "C")


GRdrawDRC(output1, points =F)

GRbox(output1, metric ='GR50', groupVariable = c('cell_line'), 
      pointColor = c("agent"))

GRbox(output1, metric ='GR50', groupVariable = c('agent'), 
      pointColor = c("cell_line"))


output2 <- GRmetrics::GRgetMetrics(output1)

cc <- scales::seq_gradient_pal("blue", "red")(seq(0,1,length.out= 7))

#cc <- grDevices::colorRampPalette(colors = c("blue", "red"))

data_figure <- subset(data_corrected, Drug %in% c("DMSO",drugs_in_screen))

data_figure <- data_figure[!(as.character(data_figure$Final_conc_uM) %in% c("33.3", "3.33")),]

data_figure$Final_conc_uM <- ifelse(data_figure$Final_conc_uM == 333, 0, data_figure$Final_conc_uM)


ggplot(subset(subset(data_figure, cell == "NCIH460"), Time > 24 & Time < 72 ) , aes(x = Time, y = Conf, color = factor(Final_conc_uM)))+
  geom_smooth()+
  facet_grid(~cell)+
  scale_color_manual(values = cc)+
  guides(color=guide_legend(title="Concentration (uM)"))


plot(output2$time, output2$GR50, xlab = "Time (h)", ylab = "GR50", col = output2$cell_line)



# # GR analysis w.r.t. time -----------------------------------------------


# heatmap GR over time

View(data_comb)

output2$drug_cell <- paste(output2$agent, output2$cell_line, sep = "_") 

output2$GR50 <- ifelse(output2$GR50 == Inf | output2$GR50 == -Inf, NA, output2$GR50)

ggplot2::ggplot(subset(output2, agent == "17-AAG") ,aes(x = time, y = drug_cell, fill = (GR50))) + geom_tile()+
  scale_fill_continuous(na.value = "white")+
  theme(axis.text=element_text(size=20), legend.text = element_text(size = 20), legend.title=element_text(size=20), axis.title = element_text(size=20))

library(drc)

### MAURO ###
### MAURO ###
### MAURO ###

#TODO create a model that we can compare

model_growth <- function(cell_tz, t, t_onset, k, maxeff, halfeff, conc, hill){
  
  maxeff <- ifelse(t_onset > t, 0, maxeff)
  
  if(t >= t_onset){
    
    cell_tz <- cell_tz * exp( ( t_onset * k ) )
    
    t <- t - t_onset
    
    return( cell_tz * exp( ( t * k *( 1 - ( (maxeff*conc ** hill )/ (halfeff + conc ** hill) ) ) ) ) )

  }else{
    
    return( cell_tz * exp( ( t * k *( 1 - ( (maxeff*conc ** hill )/ (halfeff + conc ** hill) ) ) ) ) )

    }
}
  

colfun <- colorRampPalette(c("firebrick1","firebrick4"))

time <- seq(0,3,0.01)

conc <- seq(0.5,1.5,0.2)

t_onset <- 2.5

matrix_data <- sapply(time, function(xtime) model_growth(cell_tz = 5,t = xtime, t_onset = Inf, k =  0.5, maxeff = 1, halfeff = 1.5, conc = 1.5, hill = 1.6))

matrix_data <- base::rbind(matrix_data,
                           sapply(time, function(xtime) model_growth(cell_tz = 5,t = xtime, t_onset = t_onset, k =  0.5, maxeff = 1.5, halfeff = halfeff, conc = conc, hill = 1.6)))

colnames(matrix_data) <- time

matrix_data <- t(matrix_data)

plot(rownames(matrix_data), matrix_data[,1] , type = "l", col = "blue", main = "Dynamic Growth Response", xlab = "Time [Days]", ylab = "Cell Count",
     lwd = 3)

abline(v = t_onset, col = "orange", lwd = 3)

invisible(
lapply(2:ncol(matrix_data), function (col) lines(rownames(matrix_data), matrix_data[,col], col = colfun(ncol(matrix_data))[col], lwd = 2.5))
)
lines(rownames(matrix_data), matrix_data[,1] , type = "l", col = "blue", lwd = 3)
legend(as.numeric(max(rownames(matrix_data)))/20, max(matrix_data), 
       legend = c("Control", paste("c = ", conc)), col = c("blue", colfun(ncol(matrix_data)-1)), lty = 1, lwd = 3,
       box.lty = 0)

legend(as.numeric(max(rownames(matrix_data)))/20, max(matrix_data)/1.4, 
       legend = c(paste0("D = ", t_onset)), col = c("orange"), lty = 1, lwd = 3,
       box.lty = 0)


#TODO test GR50, IC50 and GI50 for the previous data 

#TODO sample the space of growth rates from 1.9 to 3.9


#TODO create a sampling model, where we will get the best estimate and use less computer 
#TODO compare if we see any trends between drugs, which could be linked to mechanism of action. Eg. pick a drug with a very early vs. late response and compare


### UWE ###
### UWE ###
### UWE ###  

#TODO See how many cell lines change from sensitive to resistant
#TODO compare established metrics (GR50m EC50, GI50, IC50, etc.) and report to Uwe
#TODO try to prove it functionally using my data (?) on those examples
#TODO 
#TODO 

#TODO use Resazurin to prove its indeed the breaking point 

### MATTIA ###
### MATTIA ###
### MATTIA ###

#TODO  remove confluence offset by normalizing by confluence 
#TODO the time point D doesn’t match figure.. figure out a way to improve

