
library(ggplot2); library(GRmetrics)

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

data_corrected$HCT15_P1.txt <- data$HCT15_P1.txt  # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the uncorrected data

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

# plotting all data to check for inconsistencies

cc <- scales::seq_gradient_pal("lightblue", "red")(seq(0,1,length.out= length(unique(data_corrected$Final_conc_uM))))

cc <- scales::seq_gradient_pal("lightblue", "red")(seq(0,1,length.out= 9))

ggplot(subset(data_corrected, source_plate == "P1.txt" & Time <= 72 & cell == "A549" & Drug %in% c("Docetaxel", "DMSO")),
       aes(x = Time, y = Conf, color = factor(Final_conc_uM), group = factor(Final_conc_uM)))+
  geom_smooth()+
  scale_color_manual(values = cc)+
  facet_wrap(~Drug+cell)



# checking if exceptions show the same behavious as PBS or DMSO

data_corrected <- do.call(rbind, data_corrected)

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
# I will try with Clofarabine - HCT-15

# we need to select time point zero as the time we treated the samples (as per their requirements)


#time_treatment <- 24


drugs_in_screen <- c("Docetaxel")  #c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO", NA, "Water", "WATER", "H2O", "Dmso", "Control", "Ctrl"))]))

data_corrected[data_corrected$Drug %in% drugs_in_screen,]

matrix()   # columns will be drugs, rows will be cell lines

data_GRmetrics <- data_corrected

cell_line <- "HCT15"
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

data_Grmetrics_Ctr$Time <- ifelse(data_Grmetrics_Ctr$Time == min(data_Grmetrics_Ctr$Time), 0, data_Grmetrics_Ctr$Time) 
data_Grmetrics_Ttm$Time <- ifelse(data_Grmetrics_Ttm$Time == min(data_Grmetrics_Ttm$Time), 0, data_Grmetrics_Ctr$Time)

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


GRbox(output1, metric ='GR50', groupVariable = c('cell_line'), 
      pointColor = c("time"))

output2 <- GRmetrics::GRgetMetrics(output1)

ggplot(data_comb, aes(x = time, y = cell_count, color = factor(concentration)))+
  geom_point()+
  facet_wrap(~cell_line + agent)


#GI50

# Compare metrics with DTP data, and based on that choose the best metric

########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############
########## old code ##############

analyze_growth_zero <-
  function(data_source){

  data_source <- filtered_data[filtered_data$cell == cell & filtered_data$well == well,c("time", "conf")]

  # modelling confluence over time ------------------------------------------

  data <- list()

  data$store_vars <- data.frame(time = data_source[1],
                                confluence = data_source[2])

  names(data$store_vars) <- c("time", "confluence")

  # variables

  poly_degree = 5 # degrees of polynomial fit

  # modelling confluence wtih a polynomial regression model ---------

  data$growth_model <- lm(confluence ~ poly(time, degree = poly_degree, raw = t), data$store_vars)

  # we can plot it to see if the fit was good

  data$store_vars$fitted <- fitted(data$growth_model) #this function gives the predicted confluence based on the original time vector supplied

  #plot(data$store_vars$time, data$store_vars$confluence, type = "p") # plotting confluence over time
  #lines(data$store_vars$time, data$store_vars$fitted, type = "l", col = "red") #plotting fitted line

  polynomial_deriv <- function(x) {  #this function takes the first derivative, or the growth rate equation

      terms <- coef(x)

      stopifnot(names(terms)[1]=="(intercept)")

      filter_terms <- terms[-1] #remove intercept

      stopifnot(all(grepl("^poly", names(filter_terms))))

      degrees <- as.numeric(gsub(pattern = "poly\\(.*\\)", replacement = "", names(filter_terms)))

      terms_deriv <- setnames(c(filter_terms * degrees, 0), names(terms))

      terms_deriv[is.na(terms_deriv)] <- 0

      return(matrix(terms_deriv, ncol = 1))

    }

  data$store_vars$growth_rate <-
    model.matrix(data$growth_model) %*% polynomial_deriv(data$growth_model)

  return(diff(sign(data$store_vars$growth_rate)))

}

# fixing time problems in dataset ---------------------------------------

# hct15 did not record some wells, so i am correcting for this

data$hct15_p1.txt$time <- data$hct15_p1.txt$time +
  (max(data$a549_p1.txt$time) - max(data$hct15_p1.txt$time))

# filtering wells ---------------------------------------------------------

#removing bad wells

store_outlier <-

  sapply(names(data), use.names = t, function(x){

    plate_name = x

    x = data[[x]]

    good_wells <-
    filter_growth_outliers(plate_name = plate_name,data = x, time_control = 24, save_diag_plots = t,
                           save_plots_directory = "\\\\imsbnas.d.ethz.ch/sauer1/users/mauro/cell_culture_data/190310_pre_screen_test_cpddilution/figures")

    new_data <-
    x[!(x[,"well"] %in% good_wells$wells_exception),]

    return(list(new_data, exception_number = good_wells$wells_exception_number))

    })

store_outlier[[1]]

#lapply(filtered_data, function(x) plot_multi_well(input_df =  x, max_timepoint = 24))


# removing wells that echo did not transfer

source("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cell_culture_data\\190310_pre_screen_test_cpddilution\\exceptions\\log_processing.r")

source_plates <- data.frame(filenames = list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cell_culture_data\\190310_pre_screen_test_cpddilution\\20190513\\results"))

source_plates$sourceid <- ifelse(grepl(pattern = "p1", source_plates$filenames), "1msp001", "2msp001")

filtered_data <-

  lapply(names(data), function(x){

    #x = "a549_p1.txt"

    source_id <- unlist(subset(source_plates, filenames == x, sourceid))

    exception_wells <- subset(exceptions, destination.plate.name  == source_id, destination.well, drop = t)

    tmp_data <- filtered_data[[x]]

    tmp_data <- tmp_data[!(tmp_data$well %in% exception_wells),]

    tmp_data$cell <- strsplit(x, split = "_")[[1]][1]

    tmp_data$source_plate <- strsplit(x, split = "_")[[1]][2]

    return(tmp_data)

  })


# combining metadata of source.plates into growth_data

source_layout <-# import source plate layout data

  lapply(list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cpd_data\\190508_msp_layout",
                    pattern = "randomized",
                    full.names = t),
         function(x) {

           tmp_data  = readxl::read_xls(x)

           tmp_data$well = paste0(tmp_data$row, tmp_data$column)

           return(tmp_data)

         })


filtered_data <-

lapply(filtered_data, function(x){

  stopifnot(require(dplyr))

  #x = filtered_data[[1]]

  if(grepl(pattern = "p1", x = x$source_plate[1])){

    tmp_data <- inner_join( data.frame(x), source_layout[[1]][,c("well","drug" , "final_conc_um")], by = "well")

    return(tmp_data)

  }else if(grepl(pattern = "p2", x = x$source_plate[1])){

    tmp_data <- inner_join( data.frame(x), source_layout[[2]][,c("well","drug" , "final_conc_um")], by = "well")

    return(tmp_data)

  }else{

    stop("some experimental plate could not be matched to a source plate layout.")

  }

} )



# combine dataframes

filtered_data <- do.call(rbind, filtered_data)
filtered_data <- filtered_data[!(is.na(filtered_data$drug)),]


(plot <-

ggplot(subset(filtered_data, source_plate= "p1.txt"), aes(x = time, y = conf, color = factor(final_conc_um)))+
  geom_smooth(formula = y ~ poly(x, 6), method = "lm")+
  facet_wrap(~drug + cell, nrow = 5)+
  ggplot2::theme_minimal()+
  geom_vline(xintercept = 48, colour = "red"))

ggsave(plot = plot, filename = "alldoseresponseplots.png",
       path = "\\\\imsbnas.d.ethz.ch\\sauer1\\users\\mauro\\cell_culture_data\\190310_pre_screen_test_cpddilution\\figures", device = "png",
       width = 30, height = 15)



# calculating second derivative, where is first zero?


filtered_data <- subset(filtered_data, source_plate == "p1.txt" & time <= 72)

growthsolution <- filtered_data

for(cell in unique(filtered_data$cell)){

  for(well in unique(filtered_data$well)){



    print(cell)
    print(well)

    if(nrow(filtered_data[filtered_data$cell == cell & filtered_data$well == well,c("time", "conf")])== 0){
      next()
    }

    store <-

    analyze_growth_zero(

    filtered_data[filtered_data$cell == cell & filtered_data$well == well,c("time", "conf")]
    )


    store <- store[!(is.na(store))]

    store <- sum(!(store == 0))

    growthsolution[growthsolution$well == well & growthsolution$cell == cell,"conf"] <- store

  }

}

growthsolution <- growthsolution %>% group_by(well, cell) %>% slice(1)

ggplot(growthsolution, aes(cell, drug, fill = conf))+
  geom_tile()+
  scale_fill_gradient(low = "blue", high = "red")



#TODO Prepare all data and  Fit a linear regression on GI50 and report the slope beta.. if there's a general trend, it should be picked up. Make a boxplot and report the slope on the y axis for each combination

