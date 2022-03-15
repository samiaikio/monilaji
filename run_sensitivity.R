# run_sensitivity
#
# The script runs a sensitivity analysis for a selected parameter of the
# population dynamic model "function_multispecies".

# Population dynamics are simulated to reach a stable state (optional). After
# this pre-simulation, the value of the selected parameter is changed as
# specified and the populations simulated with the new value for the selected
# length of time. At the end fo the simulation, the difference and rate of
# change in population sizes are calculated relative to the stable state or
# (alternatively) the initial state of the populations. The results are
# presented as the change and rate of change in population sizes in relation to
# the values of the sensitivity parameter.
#
# Author: Sami Aikio (sami.aikio@gmail.com)
# Finnish Museum of Natural History LUOMUS
# University of Helsinki, Finland
#
# Version 2021-01-14 (sami.aikio@luke.fi)
# Natural Resources Institute Finland


# Initial settings =============================================================


# list the required packages
package_names <- c(
  "tictoc",
  "reshape2", 
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
setwd(dirname(parent.frame(2)$ofile))     # working directory to this file
#setwd("C:/Users/03180418/Documents/Monilaji_uusi3")
options(max.print=1e4)
tic()                                     # start clock


# Parameters ===================================================================

# read the function from another file
source("core_multispecies.R")     
source("functions_collected.R")     

# population parameters in an Excel file
file_name = "../Data/Kainuu_parametrit_h4.xlsx"

# simulation parameters
stochastic      <- F      # demographic stochasticity (T=on, F=off)
num_years_pre   <- NA    # pre-simulation length in years (NA = not done)
num_years_sens  <- 50     # sensitivity simulation length in years

# sensitivity parameters
sens_num        <- 7     # number of sensitivity parameter values 
#sens_name <- "density_dep_deer"; sens_values <- seq(1e3, 1e4, length.out=sens_num)
#sens_name <- "starve_wolf"; sens_values <- seq(1e0, 1e-4, length.out=sens_num)
sens_name <- "harvest_wolf"; sens_values <- seq(0, 2, length.out=sens_num)


# plotting parameters 
census_month      <- 1      # must be one in sensitivity analyses 


# loop sensitivity parameter ===================================================

if (is.na(num_years_pre)==F) {
  
  # display on console
  print(noquote("Simulation year before sensitivity analysis:"))
  
  
  # initial values from previous run -------------------------------------------
  
  # call external function for within year dynamics
  pop_sizes <- function_multispecies(
    file_name, num_years_pre, stochastic,
    NA, NA, NA)
  
  
  # extract some parameters
  num_years       <- length(unique(pop_sizes$Year))
  species_names   <- unique(pop_sizes$Species)
  pop_names       <- unique(pop_sizes$Population)
  sex_names       <- unique(pop_sizes$Sex)
  age_names       <- unique(pop_sizes$Age_Class)
  #age_names[is.na(age_names)] <- age_names[1]
  #age_names       <- age_names[is.na(age_names)==F] # remove NAs
  
  # use old initial as a template for the new
  initial      <- data.frame(read_excel(file_name, sheet="initial"))
  #initial_new      <- data.frame(read_excel(file_name, sheet="initial"))  
  initial_new      <- initial
  
  # loop for each species, population, sex and age
  for (idx_spp in 1:length(species_names)) {
    for (idx_pop in 1:length(pop_names)) {
      for (idx_sex in 1:length(sex_names)) {
        for (idx_age in 1:length(age_names)) {
          
          #  new initial value
          idx_value <- which(
            pop_sizes$Year==num_years & pop_sizes$Month==1 &
              pop_sizes$Species==species_names[idx_spp] & 
              pop_sizes$Population==pop_names[idx_pop] &
              (pop_sizes$Sex==sex_names[idx_sex] |
                 is.na(pop_sizes$Sex)==T) &
              pop_sizes$Age_Class==age_names[idx_age] 
          )
          
          end_vals <- pop_sizes$Size[idx_value]
          end_vals[is.na(end_vals)] <- T
          if (length(end_vals)==0) end_vals <- NA
          
          # row number for the new initial value
          row_idx <- which(
            initial$species==species_names[idx_spp] & 
              initial$population==pop_names[idx_pop] & 
              (initial$sex==sex_names[idx_sex] | 
                 is.na(initial$sex==sex_names[idx_sex])==T) 
          )
          
          # column number for the new initial value
          col_idx <- 3+idx_age
          if (length(row_idx) == 1) {
            initial_new[row_idx, col_idx] <- end_vals
          }
          
        }
      }
    }
  } # end loops for species etc.
  
  # display new initial sizes on console
  print(noquote("Initial population sizes for sensitivity analysis:"))
  print(initial_new)
  
  
} # end simulation for new initial values


# initiate containers
pop_sens <- list()

# display on console
print(noquote("Values of the sensitivity parameter:"))

for (sens_idx in 1:sens_num) {
  
  
  if (is.na(num_years_pre)==T) {
    initial_new <- NA
  }
  
  # sensitivity parameter value
  sens_param <- sens_values[sens_idx]
  
  # display the running sensitivity parameter values on console
  print(sens_param)
  
  # call function with new initial state
  pop_sizes <- function_multispecies(
    file_name, num_years_sens, stochastic,
    sens_name, sens_param,  initial_new)
  
  
  # extract population sizes at the census month of the last year
  idx <- pop_sizes$Month==census_month & 
    pop_sizes$Year==num_years_sens
  
  pop_end <- pop_sizes[idx, ]
  pop_end$sens_param <- sens_param
  pop_end[, c(1, 2, 8)] <- NULL  # remove year, month and time columns
  
  # store values in a list for each sensitivity parameter
  pop_sens[[sens_idx]]       <- pop_end
  
} # end of sensitivity parameter loop

# convert results into long form -----------------------------------------------

pop <- melt(pop_sens, 
            id.vars=c("Species", "Population","Sex", "Age_Class", 
                      "sens_param", "Size"))
pop <- pop[, -7]


# growth rates =================================================================

# extract initial states 
initial <- pop_sizes[pop_sizes$Mont==0, ]
initial[, c(1, 2, 8)] <- NULL


# pool sexes and age classes 
pop_total <- aggregate(Size ~ Species + Population + sens_param, 
                       data=pop, FUN=sum)

pop_init <- aggregate(Size ~ Species + Population,
                      data = initial, FUN=sum)
colnames(pop_init)[3] <- "Size_end"


# merge initial and end population sizes
pop_total <- merge(pop_total, pop_init)


# growth rates
pop_total$change <- pop_total$Size / pop_total$Size_end
pop_total$change_rate <- log(pop_total$change) / num_years_sens


# Plot population growth by sensitivity parameter value ========================

p1 <- ggplot() + 
  geom_line(data=pop_total, 
            aes(x=sens_param, y=change, color=Species), 
            size=1, linetype="solid") +
  xlab(sens_name) +
  theme_bw()

p2 <- ggplot() + 
  geom_line(data=pop_total, 
            aes(x=sens_param, y=change_rate, color=Species),
            size=1, linetype="solid") +
  xlab(sens_name) +
  theme_bw()

#png("sens.png", width=600, height=600)
print(cowplot::plot_grid(
  p1, p2, nrow=2, 
  label_size = 10, 
  label_x = 0, 
  labels = c('(a)', '(b)')))
#dev.off()


