## This code imports and clean metabolomics data

# load packages nad definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

library(R.matlab)
source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\plate_converter.R")


#import metabolomics data
setwd(path_metabolomics_in)


# # map 96wp metabolomics to 384 source -----------------------------------
#each 384 plate has four 96 source plates, to a total of 8x 96 well plates
#used for metabolomics



# #remove echo problems ---------------------------------------------------
#remove wells that we had problems with drug transfer from echo pipetting



# # correct pipetting mistakes that can be corrected ----------------------
#when transfering source plate to 96 wp, there were a few pipetting mistakes
#all pipetting mistakes were recorded and can now be corrected


#SF539	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384
#P1 on P2Q3 row D of all CLs instead of original layout due to pipetting mistake

#TODO

#IGROV1	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384
#P1 on P2Q3 row D of all CLs instead of original layout due to pipetting mistake; 

#TODO

#MDAMB231	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384 P1
#on P2Q3 row D of all CLs instead of original layout due to pipetting mistake

#TODO

#HOP62	metabolomics: hop62 p1-q1:q4, A1 is H12, front == end

#TODO

#COLO205	metabolomics: colo 205 cl3_p2_q1 row f on e and e on f

#TODO switch rows 


# # remove pipetting mistakes that cannot be corrected  -------------------


# # filter ions based on CV -----------------------------------------------



# # filter samples based on CV --------------------------------------------


#export results