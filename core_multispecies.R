# core_multispecies
#
# Function to run a simulation model of Finnish Forest Reindeer (Deer), Moose and
# an uspeciefied group "otherspecies", which are all predated on by a Wolf
# population.
#
# Population dynamics are modelled at sub-year (e.g. montly) time steps.
# Depending on the choice of parameter values, populations may be divided into
# subpopulations that are connected by monthly dispersal.
#
# Reindeer and Moose populations are structured by sex, age class and population,
# which are allowed to differ with respect to monthly survival, fecundity and
# dispersal probability. Wolf population is structured by age group (cub,
# juvenile, adult) and population. The "otherprey" group is unstructured (i.e.
# presented only by population size) and has no dispersal between
# subpopulations.
#
# The model is constructed for the Big Game Multispecies -project of the Natural
# Resources Institute Finland (Luke).
#
# Author: Sami Aikio (sami.aikio@gmail.com)
# Finnish Museum of Natural History LUOMUS
# University of Helsinki, Finland
#
# Version 2022-02-21 (sami.aikio@luke.fi)
# Natural Resources Institute Finland


# Initial settings #############################################################

# list the required packages
package_names <- c(
  "tidyverse",
  "readxl",
  "plyr",
  "dplyr",
  "reshape")

# install and attach missing packages
new_packages <- 
  package_names[!(package_names %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages) 	# install
lapply(package_names, library, character.only = TRUE)     # attach


# multispecies model ###########################################################

function_multispecies <- function(
  file_name, num_years, stochastic, 
  sens_name, sens_values, initial_new) {
  
  # code testing
  #  sens_name=NA; sens_values=NA; initial_new=NA; yr=1; i=1; j=1
  
  
  # read parameters ============================================================
  # from excel sheets 
  
  # get initial state either from excel or from the function
  if (all(is.na(initial_new))==T) {
    initial      <- data.frame(read_excel(file_name, sheet="initial"))  
  } else {
    initial <- initial_new
  }
  
  # get the rest of the parameters
  params        	<- data.frame(read_excel(file_name, sheet="params"))  
  harvest         <- data.frame(read_excel(file_name, sheet="harvest"))  

  raw_fec_deer    <- data.frame(read_excel(file_name, sheet="fec_deer"))
  raw_surv_deer   <- data.frame(read_excel(file_name, sheet="surv_deer"))
  raw_disp_deer   <- data.frame(read_excel(file_name, sheet="disp_deer"))
  connect_deer    <- data.frame(read_excel(file_name, sheet="connect_deer"))
  
  raw_fec_moose   <- data.frame(read_excel(file_name, sheet="fec_moose"))
  raw_surv_moose  <- data.frame(read_excel(file_name, sheet="surv_moose"))
  raw_disp_moose  <- data.frame(read_excel(file_name, sheet="disp_moose"))
  connect_moose   <- data.frame(read_excel(file_name, sheet="connect_moose"))
  
  raw_fec_wolf    <- data.frame(read_excel(file_name, sheet="fec_wolf"))
  raw_surv_wolf   <- data.frame(read_excel(file_name, sheet="surv_wolf")) 
  raw_disp_wolf   <- data.frame(read_excel(file_name, sheet="disp_wolf"))
  connect_wolf    <- data.frame(read_excel(file_name, sheet="connect_wolf"))
  
  
  # extract variables ==========================================================
  # species, population and age class 
  
  # change sensitivity parameter values if multiple values are given ...........
  # check if sensitivity is on harvest parameter
  harvest_sens <- sens_name %in% 
    c("harvest_deer", "harvest_moose", "harvest_wolf")
  
  # harvest parameter sensitivity
  deer_harvest_multiplier <-  
    moose_harvest_multiplier <- 
    wolf_harvest_multiplier <- 1
  if (harvest_sens==T) {  
    if (sens_name=="harvest_deer") 	deer_harvest_multiplier   <- sens_param
    if (sens_name=="harvest_moose") moose_harvest_multiplier 	<- sens_param
    if (sens_name=="harvest_wolf") 	wolf_harvest_multiplier   <- sens_param
  }
  
  # other parameter sensitivity
  if (harvest_sens==F & all(is.na(sens_values)==F)) { 
    sens_idx <- params$name == sens_name # parameter to change
    if (all(sens_idx == F)) {  # not in the "params" sheet
      stop("Invalid sensitivity parameter name!")
    }
    params$value[sens_idx] <- sens_values  # new value to the parameter
  }
  
  # empty parameter values are interpreted as infinity
  params$value[is.na(params$value)==T] <- Inf
  
  
  # extract single parameters ..................................................
  
  # Make a variable of each row in the table of scalar valued parameters
  for (p in 1:length(params$name)) {
    assign(params$name[p], params$value[p])  # CAUTION! Dynamic varible naming!
  }
  
  # names and numbers of species, populations and months
  species_names   <- unique(initial$species)
  num_species     <- length(species_names)
  pop_names       <- unique(initial$population)
  num_pops        <- length(pop_names)
  month_names     <- unique(raw_fec_deer$month)
  num_months      <- length(month_names)
  
  # convert first columns to row names in connectance matrices
  rownames(connect_deer)    <- connect_deer[, 1];   connect_deer[, 1]  <- NULL
  rownames(connect_moose)   <- connect_moose[, 1];  connect_moose[, 1] <- NULL
  rownames(connect_wolf)    <- connect_wolf[, 1];   connect_wolf[, 1]  <- NULL
  
  # initial population sizes ...................................................
  # species names must correspond to the names on the initial sizes -table
  
  # deer
  init_deer             <- initial[initial$species=="Deer", ]
  init_deer             <- init_deer[, colSums(is.na(init_deer))==0]
  deer_numeric          <- unlist(lapply(init_deer, is.numeric))
  deer_age_classes      <- colnames(init_deer)[deer_numeric]
  
  # moose
  init_moose            <- initial[initial$species=="Moose", ]
  init_moose            <- init_moose[, colSums(is.na(init_moose))==0]
  moose_numeric         <- unlist(lapply(init_moose, is.numeric))
  moose_age_classes     <- colnames(init_moose)[moose_numeric]
  
  # wolf
  init_wolf             <- initial[initial$species=="Wolf", ]
  init_wolf             <- init_wolf[, colSums(is.na(init_wolf))==0]
  wolf_numeric          <- unlist(lapply(init_wolf, is.numeric))
  wolf_age_classes      <- colnames(init_wolf)[wolf_numeric]
  
  # otherprey
  init_otherprey        <- initial[initial$species=="otherprey", ]
  init_otherprey        <- init_otherprey[, colSums(is.na(init_otherprey))==0]
  
  
  # convert initial sizes ======================================================
  # to lists where each population is a separate item
  
  # initialize lists for initial states by species and population
  initial_deer_female     <- initial_deer_male     <- 
    initial_moose_female  <- initial_moose_male    <- 
    initial_wolf          <- initial_otherprey     <- 
    list() 
  
  # loop for populations and keep only population sizes
  for (j in 1:num_pops) { 
    initial_deer_female[[j]] <- 
      subset(init_deer, sex=="f" & population==pop_names[j])[deer_numeric]
    initial_deer_male[[j]] <- 
      subset(init_deer, sex=="m" & population==pop_names[j])[deer_numeric]
    initial_moose_female[[j]] <- 
      subset(init_moose, sex=="f" & population==pop_names[j])[moose_numeric]
    initial_moose_male[[j]] <- 
      subset(init_moose, sex=="m" & population==pop_names[j])[moose_numeric]
    initial_wolf[[j]] <- 
      subset(init_wolf, population==pop_names[j])[wolf_numeric]
    initial_otherprey[[j]] <- subset(init_otherprey, 
                                     population==pop_names[j])[-c(1:2)]
  }
  
  
  # parameters as lists ========================================================
  # to lists where each population is a spearate item
  
  # initialize lists for parameters by populations 
  fec_deer    <- surv_deer_f  <- surv_deer_m  <- disp_deer_f  <- disp_deer_m <-
    fec_moose <- surv_moose_f <- surv_moose_m <- disp_moose_f <- disp_moose_m <-
    fec_wolf  <- surv_wolf    <- disp_wolf    <- list()
  
  # turn populations into list items and first column to row names
  # this does not apply to otherprey, which has no associated parameter tables
  for (j in 1:num_pops) {
    fec_deer[[j]]     <- par2list_func(raw_fec_deer, sex=NA, pop_names[j])  
    surv_deer_f[[j]]  <- par2list_func(raw_surv_deer, sex="f", pop_names[j]) 
    surv_deer_m[[j]]  <- par2list_func(raw_surv_deer, sex="m", pop_names[j]) 
    disp_deer_f[[j]]  <- par2list_func(raw_disp_deer, sex="f", pop_names[j])
    disp_deer_m[[j]]  <- par2list_func(raw_disp_deer, sex="m", pop_names[j])
    
    fec_moose[[j]]    <- par2list_func(raw_fec_moose, sex=NA, pop_names[j])  
    surv_moose_f[[j]] <- par2list_func(raw_surv_moose, sex="f", pop_names[j]) 
    surv_moose_m[[j]] <- par2list_func(raw_surv_moose, sex="m", pop_names[j]) 
    disp_moose_f[[j]] <- par2list_func(raw_disp_moose, sex="f", pop_names[j])
    disp_moose_m[[j]] <- par2list_func(raw_disp_moose, sex="m", pop_names[j])
    
    fec_wolf[[j]]     <- par2list_func(raw_fec_wolf, sex=NA, pop_names[j])  
    surv_wolf[[j]]    <- par2list_func(raw_surv_wolf, sex=NA, pop_names[j])  
    disp_wolf[[j]]    <- par2list_func(raw_disp_wolf, sex=NA, pop_names[j])  
  }
  
  # initialize containers for all population sizes
  deer_female_all     <- deer_male_all <- 
    moose_female_all  <- moose_male_all <- 
    wolf_all          <- otherprey_all <- 
    list()
  
  
  # LOOP YEARS =================================================================
  
  for (yr in 1:num_years) {  
    
    # display counter of years only when no sensitivity analysis is done
    if (is.na(sens_values)==T) cat(yr, "\n")
    
    # initialize monthly population size containers ............................
    
    # deer female and male
    monthly_deer_f <- monthly_deer_m <- 
      data.frame(matrix(byrow=T, ncol=length(deer_age_classes),  
                        rep(NA, (num_months+1)*length(deer_age_classes))))
    rownames(monthly_deer_f) <- c("initial", month_names)
    colnames(monthly_deer_f) <- deer_age_classes
    
    
    # moose female and male
    monthly_moose_f <- monthly_moose_m <- 
      data.frame(matrix(byrow=T, ncol=length(moose_age_classes), 
                        rep(NA, (num_months+1)*length(moose_age_classes))))
    rownames(monthly_moose_f) <- c("initial", month_names)
    colnames(monthly_moose_f) <- moose_age_classes
    
    # wolf 
    monthly_wolf <-       
      data.frame(matrix(byrow=T, ncol=length(wolf_age_classes), 
                        rep(NA, (num_months+1)*length(wolf_age_classes))))
    rownames(monthly_wolf) <- c("initial", month_names)
    colnames(monthly_wolf) <- wolf_age_classes
    
    # otherprey
    monthly_otherprey <-
      data.frame(matrix(byrow=T, ncol=1, rep(NA, (num_months+1))))
    rownames(monthly_otherprey) <- c("initial", month_names)
    colnames(monthly_otherprey) <- "size"
    
    # replicate container for each population 
    deer_female   <- deer_male    <- rep(list(monthly_deer_f), num_pops)
    moose_female  <- moose_male   <- rep(list(monthly_moose_f), num_pops)
    wolf                          <- rep(list(monthly_wolf), num_pops)
    otherprey                     <- rep(list(monthly_otherprey), num_pops)      
    
    # put initial size at the first row for each population
    for (j in 1:num_pops) { 
      deer_female[[j]][1, ]   <- initial_deer_female[[j]] 	# Female deer
      deer_male[[j]][1, ]     <- initial_deer_male[[j]] 	  # Male deer
      moose_female[[j]][1, ]  <- initial_moose_female[[j]] 	# moose female
      moose_male[[j]][1, ]    <- initial_moose_male[[j]] 	  # moose female
      wolf[[j]][1, ] 		      <- initial_wolf[[j]] 	        # wolf
      otherprey[[j]][1, ]     <- initial_otherprey[[j]] 	  # wolf
    }
    
    # initiate containers for dispersers
    deer_f_disp <- deer_m_disp <- moose_f_disp <- moose_m_disp <- wolf_disp <- 
      list()
    
    
    # Loop months, populations =================================================
    
    # Loop for months 
    for (i in 1:num_months) {
      
      # loop for populations
      for (j in 1:num_pops) { 
        
        # predation ============================================================
        
        # prey abundance by sex (sum over age classes)
        prey_abundance <- c(
          sum(deer_female[[j]][i, ]), sum(deer_male[[j]][i, ]),
          sum(moose_female[[j]][i, ]), sum(moose_male[[j]][i, ]), 
          sum(otherprey[[j]][i, ]))
        prey_total 	<- sum(prey_abundance)  # total prey abundance
        
        # sex-specific masses
        prey_weight <- c(mass_deer_f, mass_deer_m, mass_moose_f, mass_moose_m, 
                         mass_otherprey)
        
        # Samu's "fireplace model" of predation
        prey_encounter_probability  <- prey_abundance / prey_total
        expected_meat_per_prey  <- sum(prey_encounter_probability * prey_weight)
        required_meat 				  <- req_mass_per_wolf * sum(wolf[[j]][i, ])
        total_prey_required 	  <- required_meat / expected_meat_per_prey
        predation_probability   <- min(pred_max, total_prey_required/prey_total)
        
        prey_consumed 	    <- predation_probability * prey_abundance
        deer_f_consumed     <- prey_consumed[1]
        deer_m_consumed     <- prey_consumed[2]
        moose_f_consumed    <- prey_consumed[3]
        moose_m_consumed    <- prey_consumed[4]
        otherprey_consumed  <- prey_consumed[5]
        
        
        # DEER =================================================================
        
        # deer survival --------------------------------------------------------
        # includes predation and background mortality 
        # female
        predsurv <- predated_survivors_func(
          stochastic = stochastic, 
          pop_size = deer_female[[j]][i, ],
          prob_pred = predation_probability,
          prob_surv = surv_deer_f[[j]][i, ],
          prey_consumed = deer_f_consumed)
        predated_deer_f <- unlist(predsurv[1])  # number predated
        deer_female[[j]][i+1, ] <- unlist(predsurv[2])  # number surviving
        
        # male 
        predsurv <- predated_survivors_func(
          stochastic = stochastic, 
          pop_size = deer_male[[j]][i, ],
          prob_pred = predation_probability,
          prob_surv = surv_deer_m[[j]][i, ],
          prey_consumed = deer_f_consumed)
        predated_deer_m <- unlist(predsurv[1])  # number predated
        deer_male[[j]][i+1, ] <- unlist(predsurv[2])  # number surviving
        
        
        # deer reproduction ----------------------------------------------------
        # density dependence
        deer_all <- sum(deer_female[[j]][i+1, ] + deer_male[[j]][i+1, ])
        dd <- 1 - deer_all / (density_dep_deer + deer_all)
        dd[is.na(dd)] <- 0
        
        # number of calves
        offspring <- reproduction_func(
          stochastic,
          pop_size = deer_female[[j]][i+1, ],
          prop_male = prop_male,
          fecundity = fec_deer[[j]][i, ],
          density_dependence = dd)
        num_offspring_female <- unlist(offspring[1])
        num_offspring_male <- unlist(offspring[2])
        
        # add calves to the first age class
        deer_female[[j]][i+1, 1] <- sum(num_offspring_female) + 
          deer_female[[j]][i+1, 1]
        
        deer_male[[j]][i+1, 1] <- sum(num_offspring_male) + 
          deer_male[[j]][i+1, 1]
        
        
        # deer harvest ---------------------------------------------------------
        # harvest probability for current month
        harvest_deer <- deer_harvest_multiplier * 
          harvest[i, names(harvest)=="deer"]
        
        # deer female survival from harvest
        deer_female[[j]][i+1, ] <- 
          binomial_func(stochastic, 
                        pop_size = deer_female[[j]][i+1, ],
                        prob = 1-harvest_deer)
        
        # deer male survival from harvest
        deer_male[[j]][i+1, ] <- 
          binomial_func(stochastic, 
                        pop_size = deer_male[[j]][i+1, ],
                        prob = 1-harvest_deer)
        
        # deer dispersers ------------------------------------------------------
        # female
        deer_f_disp[[j]] <- 
          binomial_func(stochastic, 
                        pop_size = deer_female[[j]][i+1, ],
                        prob = disp_deer_f[[j]][i, ])
        
        # male
        deer_m_disp[[j]] <- 
          binomial_func(stochastic, 
                        pop_size = deer_male[[j]][i+1, ],
                        prob = disp_deer_m[[j]][i, ])
        
        # MOOSE ================================================================
        
        # moose survival -------------------------------------------------------
        # includes predation and background mortality 
        # female
        predsurv <- predated_survivors_func(
          stochastic = stochastic, 
          pop_size = moose_female[[j]][i, ],
          prob_pred = predation_probability,
          prob_surv = surv_moose_f[[j]][i, ],
          prey_consumed = moose_f_consumed)
        predated_moose_f <- unlist(predsurv[1])  # number predated
        moose_female[[j]][i+1, ] <- unlist(predsurv[2])  # number surviving
        
        # male
        predsurv <- predated_survivors_func(
          stochastic = stochastic, 
          pop_size = moose_male[[j]][i, ],
          prob_pred = predation_probability,
          prob_surv = surv_moose_m[[j]][i, ],
          prey_consumed = moose_m_consumed)  # number predated
        predated_moose_m <- unlist(predsurv[1])  # number predated
        moose_male[[j]][i+1, ] <- unlist(predsurv[2])  # number surviving
        
        # moose reproduction ---------------------------------------------------
        # density dependence
        moose_all <- sum(moose_female[[j]][i+1, ] + moose_male[[j]][i+1, ])
        dd <- 1 - moose_all / (density_dep_moose + moose_all)
        dd[is.na(dd)] <- 0
        
        # number of calves
        offspring <- reproduction_func(
          stochastic,
          pop_size = moose_female[[j]][i+1, ],
          prop_male = prop_male,
          fecundity = fec_moose[[j]][i, ],
          density_dependence = dd)
        num_offspring_female <- unlist(offspring[1])
        num_offspring_male <- unlist(offspring[2])
        
        # add calves to the first age class
        moose_female[[j]][i+1, 1] <- sum(num_offspring_female) + 
          moose_female[[j]][i+1, 1]
        
        moose_male[[j]][i+1, 1] <- sum(num_offspring_male) + 
          moose_male[[j]][i+1, 1]
        
        # moose harvest --------------------------------------------------------
        # harvest probability for current month
        harvest_moose <- moose_harvest_multiplier * 
          harvest[i, names(harvest)=="moose"]
        
        # moose female survival from harvest
        moose_female[[j]][i+1, ] <- 
          binomial_func(stochastic, 
                        pop_size = moose_female[[j]][i+1, ],
                        prob = 1-harvest_moose)
        
        # moose male survival from harvest
        moose_male[[j]][i+1, ] <- 
          binomial_func(stochastic, 
                        pop_size = moose_male[[j]][i+1, ],
                        prob = 1-harvest_moose)
        
        # moose dispersers -----------------------------------------------------
        # moose female
        moose_f_disp[[j]] <- 
          binomial_func(stochastic, 
                        pop_size = moose_female[[j]][i+1, ],
                        prob = disp_moose_f[[j]][i, ])
        
        # moose male
        moose_m_disp[[j]] <- 
          binomial_func(stochastic, 
                        pop_size = moose_male[[j]][i+1, ],
                        prob = disp_moose_m[[j]][i, ])
        
        # OTHER PREY ===========================================================
        
        # predated
        predsurv <- predated_survivors_func(
          stochastic, 
          pop_size = otherprey[[j]][i, ],
          prob_pred = predation_probability,
          prob_surv = 1, 
          prey_consumed = otherprey_consumed)
        predated_otherprey <- as.numeric(unlist(predsurv[1]))  
        otherprey[[j]][i+1, ] <- as.numeric(unlist(predsurv[2]))
        
        # otherprey renewal ....................................................
        
        # Beverton-Holt model
        if (stochastic==T) {
          otherprey[[j]][i+1, ] <- sum(
            rpois(n = otherprey[[j]][i+1, ], 
                  lambda = r_otherprey) * 
              rbinom(
                n = otherprey[[j]][i+1, ], 
                size = 1, 
                prob = 1/(1+otherprey[[j]][i+1, ] / 
                            (K_otherprey/(r_otherprey-1)))))
        } else {
          otherprey[[j]][i+1, ] <- r_otherprey * otherprey[[j]][i+1, ] /
            (1 + otherprey[[j]][i+1, ] / (K_otherprey/(r_otherprey-1)))
        }
        otherprey[[j]][i+1, is.na(otherprey[[j]][i+1, ])] <- 0
        
        
        # WOLF  ================================================================
        
        # total prey caught by wolf
        num_predated <- sum(predated_deer_f + predated_deer_m) + 
          sum(predated_moose_f + predated_moose_m) + 
          as.numeric(predated_otherprey)
        
        mass_predated <- 
          sum(mass_deer_f * predated_deer_f + mass_deer_m * predated_deer_m) + 
          sum(mass_moose_f * predated_moose_f + 
                mass_moose_m * predated_moose_m) + 
          mass_otherprey * as.numeric(predated_otherprey)
        
        # wolf starvation due to lack of prey
        surv_starvation <- exp(-starve_wolf * (wolf[[j]][i, ] / mass_predated))
        #surv_starvation <- min(mass_predated / required_meat, 1)
        surv_starvation = 1  # no starvation
        
        surv_starvation[is.na(surv_starvation)] <- 0
        
        # wolf survival --------------------------------------------------------
        wolf[[j]][i+1, ] <- binomial_func(
          stochastic,
          pop_size = wolf[[j]][i, ],
          prob = surv_starvation * surv_wolf[[j]][i, ] )
        
        # wolf reproduction ----------------------------------------------------
        offspring <- reproduction_func(
          stochastic,
          pop_size = wolf[[j]][i+1, ],
          prop_male = prop_male,
          fecundity = fec_wolf[[j]][i, ],
          density_dependence = dd)
        num_offspring_female <- unlist(offspring[1])
        num_offspring_male <- unlist(offspring[2])
        
        wolf[[j]][i+1, 1] <- wolf[[j]][i+1, 1] + 
          sum(num_offspring_female + num_offspring_male)
        
        # wolf maturation ------------------------------------------------------
        # number of juveniles maturing to reproductive stage 
        if (stochastic==T) {
          num_mature <- sum(rbinom(
            n = 1,
            size = wolf[[j]][i+1, 2],  # juveniles 
            prob = wolf_maturation))
        } else {
          num_mature <- wolf_maturation * wolf[[j]][i+1, 2]
        }
        
        # vacant territories is a ceiling for maturing wolf
        vacant_territories <- 
          max(0, wolf_territories - as.numeric(tail(wolf[[j]][i+1, 3])))
        num_mature <- min(c(num_mature, vacant_territories))
        
        # update population sizes after maturation
        wolf[[j]][i+1, 2] <- wolf[[j]][i+1, 2] - num_mature
        wolf[[j]][i+1, 3] <- wolf[[j]][i+1, 3] + num_mature 
        
        # wolf harvest ---------------------------------------------------------
        
        # harvest probability for current month
        harvest_wolf <- wolf_harvest_multiplier * 
          harvest[i, names(harvest)=="wolf"]
        
        # wolf survival from harvest
        wolf[[j]][i+1, ] <- binomial_func(
          stochastic,
          pop_size = wolf[[j]][i+1, ],
          prob = 1-harvest_wolf)
        
        # wolf dispersers ------------------------------------------------------
        wolf_disp[[j]] <- binomial_func(
          stochastic,
          pop_size = wolf[[j]][i+1, ],
          prob = disp_wolf[[j]][i, ])
        
        # remove dispersers ====================================================
        # from the source population 
        if (num_pops > 1) { # skip dispersal if only one population
          deer_female[[j]][i+1, ] <- deer_female[[j]][i+1, ] - deer_f_disp[[j]]
          deer_male[[j]][i+1, ]   <- deer_male[[j]][i+1, ] - deer_m_disp[[j]]
          moose_female[[j]][i+1, ]<- moose_female[[j]][i+1, ] -moose_f_disp[[j]]
          moose_male[[j]][i+1, ]  <- moose_male[[j]][i+1, ] - moose_m_disp[[j]]
          wolf[[j]][i+1, ]        <- wolf[[j]][i+1, ] - wolf_disp[[j]]
        }
        
      } # end loop for populations .............................................
      
      
      # distribute dispersers ==================================================
      # to sink populations
      if (num_pops > 1) {  # skip dispersal if only one population
        for (source_idx in 1:num_pops) {  # source population
          for (sink_idx in 1:num_pops) {  # target population
            
            deer_female[[sink_idx]][i+1, ] <- deer_female[[sink_idx]][i+1, ] + 
              connect_deer[source_idx, sink_idx] * deer_f_disp[[source_idx]]
            
            deer_male[[sink_idx]][i+1, ] <- deer_male[[sink_idx]][i+1, ] + 
              connect_deer[source_idx, sink_idx] * deer_m_disp[[source_idx]]
            
            moose_female[[sink_idx]][i+1, ] <- moose_female[[sink_idx]][i+1, ] + 
              connect_moose[source_idx, sink_idx] * moose_f_disp[[source_idx]]
            
            moose_male[[sink_idx]][i+1, ] <- moose_male[[sink_idx]][i+1, ] + 
              connect_moose[source_idx, sink_idx] * moose_m_disp[[source_idx]]
            
            wolf[[sink_idx]][i+1, ] <- wolf[[sink_idx]][i+1, ] + 
              connect_wolf[source_idx, sink_idx] * wolf_disp[[source_idx]] 
            
          }
        }
      }
      
    } # end of loop for months .................................................
    
    # end of a year management =================================================
    # book keeping to move populations to the next year
    
    # everything done one population at the time
    for (j in 1:num_pops) {
      
      # Add month, year & population names in the results ......................
      
      # remove initial sizes from populations from all but the first year
      if (yr > 1) {
        deer_female[[j]]  <- tail(deer_female[[j]], -1)
        deer_male[[j]]    <- tail(deer_male[[j]], -1)
        moose_female[[j]] <- tail(moose_female[[j]], -1)
        moose_male[[j]]   <- tail(moose_male[[j]], -1)
        wolf[[j]]         <- tail(wolf[[j]], -1)
        otherprey[[j]]    <- tail(otherprey[[j]], -1)
        first_month       <- 1 
      } else {
        first_month       <- 0 # initial state denoted as Month=0
      }
      
      # add columns for month numbers, years and populations ...................
      deer_female[[j]]$Month          <- deer_male[[j]]$Month   <- 
        moose_female[[j]]$Month       <- moose_male[[j]]$Month  <- 
        wolf[[j]]$Month               <- otherprey[[j]]$Month   <- 
        first_month:num_months
      
      rep_months = length(first_month:num_months) # number of repetitions 

      deer_female[[j]]$Year           <- deer_male[[j]]$Year   <- 
        moose_female[[j]]$Year        <- moose_male[[j]]$Year  <- 
        wolf[[j]]$Year                <- otherprey[[j]]$Year   <- 
        rep(yr, rep_months) 
      deer_female[[j]]$Population     <- deer_male[[j]]$Population   <- 
        moose_female[[j]]$Population  <- moose_male[[j]]$Population  <- 
        wolf[[j]]$Population          <- otherprey[[j]]$Population   <- 
        rep(pop_names[j], rep_months)
      
      # promote to next age class ----------------------------------------------
      
      # deer and moose
      initial_deer_female[[j]]  <- promote_func(pop_size = deer_female[[j]])
      initial_deer_male[[j]]    <- promote_func(pop_size = deer_male[[j]])
      initial_moose_female[[j]] <- promote_func(pop_size = moose_female[[j]])
      initial_moose_male[[j]]   <- promote_func(pop_size = moose_male[[j]])
      
      # wolf 
      if (length(wolf_age_classes) > 1) {
        juveniles <- sum(tail(wolf[[j]], 1)[1:2]) # promoted plus surviving
        pop_oldest 	<- as.numeric(tail(wolf[[j]], 1)[length(wolf_age_classes)])
        initial_wolf[[j]] <- c(0, juveniles, pop_oldest) # first age class zero
      } else {
        initial_wolf[[j]] <- tail(wolf[[j]], 1)
      }
      
      # otherprey 
      # (no age classes, i.e., only one value for population size)
      initial_otherprey[[j]] <- tail(otherprey[[j]], 1)[, 1]
      
    } # end loop for populations 
    
    # add results to containers (list of lists): [[year]][[population]]
    deer_female_all[[yr]]   <- deer_female
    deer_male_all[[yr]]     <- deer_male
    moose_female_all[[yr]]  <- moose_female
    moose_male_all[[yr]]    <- moose_male
    wolf_all[[yr]]          <- wolf
    otherprey_all[[yr]]     <- otherprey
    
  } # end years
  
  
  # reshape in long format =====================================================
  # (each variable in its own column) 
  
  # deer 
  deer_female_long <- longformat_func(pop_size = deer_female_all, sex = "F")
  deer_male_long <- longformat_func(pop_size = deer_male_all, sex = "M")
  
  # moose
  moose_female_long <- longformat_func(pop_size = moose_female_all, sex = "F")
  moose_male_long <- longformat_func(pop_size = moose_male_all, sex = "M")
  
  # wolf 
  wolf_long <- longformat_func(pop_size = wolf_all, sex = "NA")
  wolf_long$Species <- "Wolf"
  
  # otherprey 
  otherprey_long <- longformat_func(pop_size = otherprey_all, sex = "NA")
  otherprey_long$Age_Class <- deer_age_classes[1]
  otherprey_long$Species <- "otherprey"
  
  # output ---------------------------------------------------------------------  
  
  # combine sexes in the same dataframe 
  deer_female_long$Sex  <- "f"
  deer_male_long$Sex    <- "m"
  deer_long <- rbind(deer_female_long, deer_male_long)
  deer_long$Species <- "Deer"
  
  moose_female_long$Sex   <- "f"
  moose_male_long$Sex     <- "m"
  moose_long <- rbind(moose_female_long, moose_male_long)
  moose_long$Species <- "Moose"
  
  # make sure all species have the same columns in the same order 
  moose_long <- moose_long[, names(deer_long)]
  wolf_long <- wolf_long[, names(deer_long)]
  otherprey_long <- otherprey_long[, names(deer_long)]
  
  # combine results in the same dataframe
  pop <- rbind(deer_long, moose_long, wolf_long, otherprey_long)
  
  # linear time axis for plotting 
  pop$Time      <- pop$Year - 1 + pop$Month/num_months
  
  # output  
  return(pop)
  
} # end of function

