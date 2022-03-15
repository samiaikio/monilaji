# functions_collected.R ########################################################
#
# This is a collection of functions called inside function_multispecies.
#
# sami.aikio@luke.fi (2022)


# predated_survivors_func ######################################################
# Returns the numbers of prey animals that are
# 1. predated
# 2. surviving both predation and background mortality

predated_survivors_func <- function(
  stochastic, 
  pop_size, 
  prob_pred, 
  prob_surv,
  prey_consumed) {
  
  # make a data frame of scalar population size (just to give it a column name)
  if (is.list(pop_size)==F) pop_size <- data.frame(N=pop_size)
  
  age_classes     <- names(pop_size)
  num_age_classes <- length(age_classes)
  pop_size        <- as.numeric(pop_size)
  prob_surv       <- as.numeric(prob_surv)
  
  if (stochastic==T) {  # stochastic predation
    # numbers predated
    predated_num <- 
      data.frame(rbind(rbinom(num_age_classes, pop_size, prob_pred)))
    colnames(predated_num) <- age_classes
    
    # numbers surviving background and predation mortality
    survived_num <- 
      data.frame(rbind(rbinom(num_age_classes, 
                              pop_size - as.numeric(predated_num), 
                              prob_surv)))
    colnames(survived_num) <- age_classes
    
  } else {  # deterministic predation
    predated_num <- prey_consumed * pop_size / sum(pop_size)
    survived_num <- prob_surv * (pop_size - predated_num)
    
  }
  
  return(list(predated_num, survived_num))
}


# reproduction_func ############################################################

reproduction_func <- function(
  stochastic,
  pop_size,
  prop_male,
  fecundity,
  density_dependence) {
  
  num_age_classes <- length(names(pop_size))
  pop_size        <- as.numeric(pop_size)
  fecundity       <- as.numeric(fecundity)
  
  # number of offspring produced in each each age class
  num_offspring_female <- num_offspring_male <- NULL
  
  if (stochastic==T) {
    for (idx in 1:num_age_classes) {
      # reproduction only if age class is not zero
      if (pop_size[idx] > 0) {
        num_offspring <- sum(rpois(
          n = as.numeric(pop_size[idx]),
          lambda = fecundity[idx] * density_dependence))
        num_offspring_male[idx] <- rbinom(1, num_offspring, prop_male)
        num_offspring_female[idx] <- num_offspring - num_offspring_male[idx]
      }  
    } 
  } else {  # deterministic
    num_offspring_female  <- (1-prop_male) * sum(fecundity * pop_size) *
      density_dependence
    num_offspring_male  <- prop_male * sum(fecundity * pop_size) *
      density_dependence
  }
  num_offspring_female[is.na(num_offspring_female)] <- 0
  num_offspring_male[is.na(num_offspring_male)] <- 0
  
  return(list(num_offspring_female, num_offspring_male))
}


# binomial_func ################################################################

binomial_func <- function(
  stochastic,
  pop_size,
  prob) {
  
  age_classes     <- names(pop_size)
  num_age_classes <- length(age_classes)
  pop_size        <- as.numeric(pop_size)
  prob            <- as.numeric(prob)
  
  if (stochastic==T) {
    occurrences <- 
      data.frame(rbind(rbinom(num_age_classes, pop_size, prob)))
    colnames(occurrences) = age_classes
  } else {
    occurrences <- prob * pop_size
  }
  occurrences[is.na(occurrences)] <- 0
  return(occurrences)
}


# promote_func #################################################################
# at the end of a year, promote individuals to the next age class

promote_func <- function(pop_size) {
  
  num_age_classes <- sum(!(names(pop_size) %in% 
                             c("Month", "Year", "Population")))
  
  if (num_age_classes > 1) {
    pop_promote <- as.numeric(tail(pop_size, 1)[1:num_age_classes-1])
    pop_oldest 	<- as.numeric(tail(pop_size, 1)[num_age_classes])
    
    pop_promote[num_age_classes-1] <- 
      pop_promote[num_age_classes-1] + pop_oldest # keep survivors
    pop_initial <- c(0, pop_promote)  # first age class at zero
  } else {
    pop_initial <- as.numeric(tail(pop_size, 1))[1]
  }
  pop_initial[is.na(pop_initial)] <- 0
  return(pop_initial) 
}


# longformat_func ##############################################################
# reshape values in long format (one variable per column)

longformat_func <- function(pop_size, sex) {
  
  pop_size_long <- melt(pop_size, id.vars=c("Year", "Month", "Population"))
  pop_size_long <- pop_size_long[, -c(6, 7)]
  colnames(pop_size_long)  <- 
    c("Year", "Month", "Population", "Age_Class", "Size")
  pop_size_long$Sex = sex
  return(pop_size_long)
}


# par2list_func ################################################################

par2list_func <- function(
  raw_params, 
  sex,
  pop_names) {
  
  if (is.na(sex)) {  # sex independent parameters
    neat_params <- raw_params[raw_params$population == pop_names, ]
    rownames(neat_params)      <- neat_params[, 1]
    neat_params[, c(1:2)]      <- NULL
  } else {  # sex dependent parameters
    neat_params <- raw_params[
      raw_params$sex == "f" & raw_params$population == pop_names, ]
    rownames(neat_params)   <- neat_params[, 1]
    neat_params[, c(1:3)]   <- NULL
  }
  return(neat_params)
}





