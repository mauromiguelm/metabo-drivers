
library(ggplot2); library(RColorBrewer)

# import functions

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\CellCultureAnalysis.R")

source("C:\\Users\\masierom\\polybox\\Programing\\Project_exometabolites\\modelling_growth_curves.R")

# Importing exceptions

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\exceptions\\log_processing.r")

## Importing source plates

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData")

source_plates <- data.frame(
  filenames = list.files(
    path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\Growthdata\\",
    pattern = "[P][1-2][.txt]",
    recursive = T))

# source_plates <- rbind(
#   data.frame(filenames = list.files ("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData\\20190828\\results", 
#                                                    pattern = ".txt",
#                                                    recursive = T)))

source_plates$sourceid <- ifelse(grepl(pattern = "P1", source_plates$filenames), "1MSP001", "2MSP001")

source_plates$batch <- ifelse(grepl(pattern = "20190513", source_plates$filenames), "batch_1", "batch_2")

#import growth curves

fileNames <- as.character(source_plates$filenames)

data <- 
  
  lapply(fileNames, function(x) {
    
    #x = fileNames[1]
    
    df <-
      ReadIncuCyteData(read_platemap = F, FileName_IncuCyte = x, Plate_size = 384, FileDirectory = getwd())
    
    df$Time <- as.numeric(df$Time)
    
    return(df)
    
  })

names(data) <- fileNames


#FIXME ACHN and M14 should have time - 24h, and then exclude the negative times

data$`20190828/results/ACHN_CL3_P1.txt`$Time <- data$`20190828/results/ACHN_CL3_P1.txt`$Time -24

data$`20190828/results/ACHN_CL3_P1.txt` <- subset(data$`20190828/results/ACHN_CL3_P1.txt`, Time >= 0)

data$`20190828/results/ACHN_CL3_P2.txt`$Time <- data$`20190828/results/ACHN_CL3_P2.txt`$Time -24

data$`20190828/results/ACHN_CL3_P2.txt` <- subset(data$`20190828/results/ACHN_CL3_P2.txt`, Time >= 0)

data$`20190828/results/M14_CL2_P1.txt`$Time <- data$`20190828/results/M14_CL2_P1.txt`$Time -24

data$`20190828/results/M14_CL2_P1.txt` <- subset(data$`20190828/results/M14_CL2_P1.txt`, Time >= 0)

data$`20190828/results/M14_CL2_P2.txt`$Time <- data$`20190828/results/M14_CL2_P2.txt`$Time -24

data$`20190828/results/M14_CL2_P2.txt` <- subset(data$`20190828/results/M14_CL2_P2.txt`, Time >= 0)


####

data$`20191029/results/EKVX_CL2_P1.txt`$Time <- data$`20191029/results/EKVX_CL2_P1.txt`$Time -24

data$`20191029/results/EKVX_CL2_P1.txt` <- subset(data$`20191029/results/EKVX_CL2_P1.txt`, Time >= 0)

data$`20191029/results/EKVX_CL2_P2.txt`$Time <- data$`20191029/results/EKVX_CL2_P2.txt`$Time -24

data$`20191029/results/EKVX_CL2_P2.txt` <- subset(data$`20191029/results/EKVX_CL2_P2.txt`, Time >= 0)



# correcting time for HCT-15: HCT-15 did not record some wells, so I am correcting for this

data_corrected <- data

data_corrected$`20190513/results/HCT15_P1.txt`$Time <- data_corrected$`20190513/results/HCT15_P1.txt`$Time +
    (max(data_corrected$`20190513/results/A549_P1.txt`$Time) - max(data_corrected$`20190513/results/HCT15_P1.txt`$Time))


# removing bad wells based on first 24h  growth

skip_outlier <- c("20190513/results/HCT15_P1.txt", # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data
                  "20190828/results/ACHN_CL3_P1.txt", # only 24h of growth
                  "20190828/results/ACHN_CL3_P2.txt") # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data
 

data_skip <- data_corrected[skip_outlier]

data_corrected <- data_corrected[!(names(data_corrected) %in% skip_outlier)]

tmp_names <- names(data_corrected)

data_corrected <-
  
  lapply(names(data_corrected), FUN = function(x){
    
    print(x)
    
    # plate_idx = names(data_corrected)[14]   #FIXME the name contains / which is not allowed as a finame.. use regex to only get the real plate name or remove the backslashed
    
    # x = names(data_corrected)[14]
    
    plate_name = strsplit(x, split = "/")[[1]][3]
    
    x = data_corrected[[x]]
    
    good_wells <-
      filter_growth_outliers(plate_name = plate_name, data = x, time_control = 24, save_diag_plots = T,
                             save_plots_directory = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/Mauro/cell_culture_data/190310_LargeScreen/figures")
    
    new_data <-
      x[!(x[,"Well"] %in% good_wells$wells_exception),]
    
    print(paste("Plate:", plate_name,
                "Exception-Number:",good_wells$wells_exception_number
    ))
    
    return(new_data)
    
  })


names(data_corrected) <- tmp_names

data_corrected <- append(data_corrected, data_skip)

data_corrected <- data_corrected[fileNames]

rm(data_skip, tmp_names)

# Fit the confluence for each well, and return fitted confluence. Keep the same plate map structure.

base::lapply(names(data_corrected), function(plate_name){
  
  #plate_name =  names(data_corrected)[3] #FIXME delete me
  
  r.2.threshold = 0.8
  
  data_raw <- data_corrected[[plate_name]]
  
  plate_name <- strsplit(plate_name, split = "/")[[1]][3]
  
  fitted_data <- 
    
    lapply(unique(data_raw$Well), function(idx_well){
      
      #idx_well <- unique(data_raw$Well)[2] #FIXME delete me 
      
      data_well_raw <- subset(data_raw[,c("Time", "Conf", "Well")], Well == idx_well)
      
      min_scan_Time <- ceiling(min(data_well_raw$Time))
      
      max_scan_Time <- floor(max(data_well_raw$Time))
      
      model_pred_poly <- get_growthMetrics(data_well_raw, degree = 7)[[3]]
      
      time_sequence <- seq(min_scan_Time ,max_scan_Time,1)
      
      model_well_pred <- predict(model_pred_poly, data.frame(Time = time_sequence))
      
      pred_conf <- data.frame(Well = idx_well, Time = time_sequence, Conf = model_well_pred)
      
      well_metrics <- summary(model_pred_poly)
      
      well_metrics <- data.frame(Well = idx_well, adj.r.2 = well_metrics$adj.r.squared)
      
      return(list(pred_conf, well_metrics))
      
    })
  
  # save diagnostic plots
  
  diagnostics_well <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[2]]) ) #for each plate, plot diagnostic plots with the R-squared.. check if any of fits have ploblems
  
  save_plots_directory = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/Mauro/cell_culture_data/190310_LargeScreen/figures"
  
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
  
  lapply(list.files("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\cpd_data\\large_screen_plate_layout",
                    pattern = "randomized",
                    full.names = T,
                    recursive = T),
         function(x) {
           
           #x =  "\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\cpd_data\\large_screen_plate_layout/190806_NewMSP_Layout/randomized_layout_1MSP_batch2.xls"
           
           tmp_data  = readxl::read_xls(x)
           
           tmp_data$Well = paste0(tmp_data$Row, tmp_data$Column)
           
           return(tmp_data)
           
         })


data_corrected <-
  
  lapply(names(data_corrected), function(x){
    
    stopifnot(require(dplyr))
    
    #x = "20190513/results/A549_P1.txt" 
    
    tmp_data = data_corrected[[x]]
    
    match_source <- source_plates[source_plates$filenames == x,c("sourceid", "batch")] 
    
    if(match_source$batch == "batch_1"){
      
      if(grepl(pattern = "1MSP", x = match_source$sourceid)){
        
        tmp_data <- inner_join( tmp_data, source_layout[[1]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else if(grepl(pattern = "2MSP", x =match_source$sourceid)){
        
        tmp_data <- inner_join(  tmp_data, source_layout[[2]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else{
        
        stop("some experimental plate could not be matched to a source plate layout.")
        
      }
      
    } else if(match_source$batch == "batch_2"){
      
      if(grepl(pattern = "1MSP", x = match_source$sourceid)){
        
      tmp_data <- inner_join( tmp_data, source_layout[[3]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
      return(tmp_data)
        
      }else if(grepl(pattern = "2MSP", x =match_source$sourceid)){
        
        tmp_data <- inner_join(  tmp_data, source_layout[[4]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else{
        
        stop("some experimental plate could not be matched to a source plate layout.")
        
      }
      
    }else{
      
      stop("some information is wrong with the match_source$source_id")
      
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
    
    tmp_data$cell<- strsplit(x, split = "/")[[1]][3]
    
    tmp_data$cell <- strsplit(tmp_data$cell, split = "_")[[1]][1]
    
    tmp_data$source_plate <- strsplit(x, split = "_")[[1]][2]
    
    return(tmp_data)
    
  })

names(data_corrected) <- fileNames

# checking if exceptions show the same behavious as PBS or DMSO

data_corrected <- do.call(rbind, data_corrected)

data_corrected <- data_corrected[!(is.na(data_corrected$Final_conc_uM)), ]

# removing DMSO that is not 33.3 and 3.33, as it causes many data analysis problems

tmp_data_corrected <- data_corrected

data_corrected <- subset(data_corrected, !(Drug == "DMSO" & !(Final_conc_uM %in% c(333,367)))) #FIXME include DMSO that increased volume

rm(exceptions, source_layout, source_plates, tmp_data_corrected, fileNames, skip_outlier)
