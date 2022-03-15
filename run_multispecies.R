# run_multispecies
#
# The script runs a population dynamic model by defining the arguments and
# calling "core_multispecies", calculates the sum of age class and sex
# specific population sizes for each species and presents the results as a time
# series.
#
# Sami Aikio (sami.aikio@luke.fi)

# Initial settings =============================================================

# list the required packages
package_names <- c(
  "tictoc",
  "ggplot2",
  "cowplot")

# install and attach missing packages
new_packages <- 
  package_names[!(package_names %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages) 	# install
lapply(package_names, library, character.only = TRUE)     # attach


# clean up and settings
if (.Device != "null device") dev.off()   # clear graphics, if any
rm(list = ls())			                      # clear workspace memory
cat("\014")  				                      # clear console
dirn <- dirname(parent.frame(2)$ofile)

setwd(dirname(parent.frame(2)$ofile))     # working directory to this file

tic()                                     # start clock

options(max.print=1000)                   # number of rows in console
options(warn=1) # display warnings (warn=2 turns warnings into errors) 



# read external files ==========================================================
source("core_multispecies.R")     
source("functions_collected.R")     


# Parameters ===================================================================

# name of the Excel file containing population parameters
file_name = "../Data/Kainuu_parametrit_h4.xlsx"

# simulation parameters
stochastic      <- T      # demographic stochasticity (T=on, F=off)
num_years       <- 50   # simulation length in years, was 20

# plotting parameters
census_month    <- 6    # month number for plotting (NA = all months plotted)

# Luke colours 
cols <- c(
  "#78be20",  # green
  "#ff8200",  # orange
  "#0033a0",  # navy
  "#54585a",  # dark gray
  "#7f3f98",  # violet
  "#e13c98",  # pink
  "#00B5E2"   # light blue
) 

# call function_multispecies ===================================================
# the actual simulation is done inside this function

# call the function
pop_sizes <- function_multispecies(
  file_name, num_years, stochastic,
  sens_name=NA, sens_values=NA, initial_new=NA)

#sens_name=NA; sens_values=NA; initial_new=NA
#pop_sizes <- function_multispecies()


# sum over age classes and sexes ===============================================

pop <- aggregate(data=pop_sizes, FUN=sum, 
                 Size ~ Population + Time + Month + Species)
head(pop)

# keep population sizes only from the annual census month
if (is.na(census_month)==F) {
  pop <- pop[pop$Month==census_month, ]
}



# Plot all species in the same graph ===========================================
# (summed over sex and age class)

ggplot(data=pop, aes(x=Time, y=Size, group=Species)) +
  geom_line(size=1, aes(color=Species)) +
  theme_bw() + 
  scale_color_manual(values=cols) +
  labs(x = "Time", y = "Population size") +
  theme(legend.position="right")



# Plot all data of one species =================================================
# (messy, if there are many age classes and populations)

if (F) {
  pop_spp <- pop_sizes[pop_sizes$Species=="Deer", ]
  ggplot(data=pop_spp, aes(x=Time, y=Size, color=Age_Class, linetype=Population)) +    
    geom_line(size=1) +
    theme_bw()  
}


# Output =======================================================================

# write results in a file 
write.csv2(pop, file="../Output/pop.csv")


toc() # stop clock

print("Done!")


