get_category_structure <- function(categoryStructure, #either Information integration or unidimensional
                                   stim_per_category=200){
  
  if(categoryStructure=="II")
    return(generate_ii_coords(stim_per_category))
  
  if(categoryStructure=="UD")
    return(generate_ud_coords(stim_per_category))
  
}

generate_ud_coords <- function(stim_per_category=200){
  
  ## generates 2D stimuli using normal distributions according to
  ## the category structure defined by define_UD_structure
  
  # require(mvtnorm)
  
  structure <- data.frame(
    "stimulusNo" = 1:(2*stim_per_category),
    "Xvalue"=NA, "Yvalue"=NA,
    "category" =  factor(rep(c(0,1), each=stim_per_category),
                         labels=c("A","B"))
  )
  
  parameters <- define_UD_structure()
  
  structure[structure$category=="A", c("Xvalue","Yvalue")] <-
    rmvnorm(
      n=stim_per_category,
      mean=c(parameters$meanA_x, parameters$meanA_y),
      sigma=diag(c(parameters$var_x, parameters$var_y))
    )
  
  structure[structure$category=="B", c("Xvalue","Yvalue")] <-
    rmvnorm(
      n=stim_per_category, 
      mean=c(parameters$meanB_x, parameters$meanB_y),
      sigma=diag(c(parameters$var_x, parameters$var_y))
    )
  
  # detach("package:mvtnorm", unload=TRUE)
  
  return(structure)
}

generate_ii_coords <- function(stim_per_category=200){
  
  ## generates 2D stimuli using normal distributions according to
  ## the category structure defined by define_UD_structure
  
  # require(mvtnorm)
  
  structure <- data.frame(
    "stimulusNo"=1:(2*stim_per_category),
    "Xvalue"=NA, "Yvalue"=NA,
    "category" = factor(rep(c(0,1), each=stim_per_category), 
                        labels=c("A","B"))
  )
  
  parameters <- define_II_structure()
  
  structure[structure$category=="A", c("Xvalue","Yvalue")] <-
    rmvnorm(
      n=stim_per_category, 
      mean=c(parameters$meanA_x, parameters$meanA_y),
      sigma=matrix(c(parameters$var_x, parameters$var_xy,
                     parameters$var_xy, parameters$var_y), 2, 2)
    )
  
  structure[structure$category=="B", c("Xvalue","Yvalue")] <-
    rmvnorm(
      n=stim_per_category, 
      mean=c(parameters$meanB_x, parameters$meanB_y),
      sigma=matrix(c(parameters$var_x, parameters$var_xy,
                     parameters$var_xy, parameters$var_y), 2, 2)
    )
  
  # detach("package:mvtnorm", unload=TRUE)
  
  return(structure)
}

define_UD_structure <- function(){
  data.frame(meanA_x=35.86, meanA_y=50, meanB_x=64.14, meanB_y=50,
             var_x=16.33, var_y=355.55, var_xy=0, angle=0)
}

define_II_structure <- function(){
  data.frame(meanA_x=40, meanA_y=60, meanB_x=60, meanB_y=40,
             var_x=185.94, var_y=185.94, var_xy=169.61, angle=pi/4)
}

## Get different responses

## Responses using a unidimensional strategy
get_ud_responses <- function(coord_vector, boundary=50,
                             perceptual_noise, boundary_noise){
  
  ## generates unidimensional responses (based on dim1)
  ## args:
  ##   coord_vector: vector of either x- or y-values (depending on
  ##   the ud version)
  ##   boundary: position of boundary
  ##   perceptual_noise: perceptual noise
  ##   boundary_noise: variation in implementing boundary (norm
  ##   distributed)
  ##   exp_noise: uniformly distributed noise respresenting
  ##   participant idiocy
  ## returns:
  ##   factor of A-B responses equal to the length of the category
  ##   input
  
  ## Setup
  
  perceived_coords <- apply_perceptual_noise(coord_vector,
                                             noise = perceptual_noise)
  
  responses <- apply_boundary(perceived_coords, boundary,
                              boundary_noise)
  
  return(factor(responses, labels=c("A","B")))
  
}

# Responses using an information integration strategy
get_ii_responses <- function(category_structure, gradient,
                             y_intercept, perceptual_noise,
                             boundary_noise){
  
  ## generates information-integration responses (based on dim1)
  ## args:
  ##   category_structure: data frame describing the category
  ##   structure
  ##   perceptual_noise: perceptual noise (assumes no
  ##     correlation between dim 1 and 2)
  ##   boundary: position of boundary
  ##   boundary_noise: variation in implementing boundary (norm
  ##   distributed)
  ##   exp_noise: uniformly distributed noise respresenting
  ##   participant idiocy
  ## returns:
  ##   vector of A-B responses equal to the length of the category
  ##   input
  
  distances <- point2line_distance(gradient, y_intercept,
                                   stim_coords=category_structure[
                                     ,c("Xvalue","Yvalue")])
  
  perceived_dist <- apply_perceptual_noise(distances,
                                           noise=perceptual_noise)
  
  
  responses <- factor(apply_boundary(perceived_dist, boundary=0,
                                     boundary_noise), 
                      labels=c("A","B"))
  
  return(correct_responses(responses, category_structure$category))
  
}

get_random_responses <- function(bias=0.5, no_trials){
  ## generates random responses
  ## args:
  ##   bias: bias towards responding A or B, 0=all A; 1=all B
  ## returns:
  ##   vector of A-B responses equal to the length of the category input
  
  if (bias>0 && bias<1){
    
    responses <- factor(sample(c(0,1), size=no_trials, replace=T,
                               prob=c(1-bias, bias)),
                        labels=c("A","B"))
    
  } else if(bias==0){
    
    responses <- factor(rep("A", no_trials))
    levels(responses)<-c(levels(responses), "B")
    
  } else if (bias==1){
    
    responses <- factor(rep("B", no_trials))
    levels(responses)<-c("A", levels(responses))
    
  }
  
  return(responses)
}

get_response_parameters <- function(simulation_frame){
  
  return(list(perceptual_noise=simulation_frame$perceptualNoise, 
              boundary_noise=simulation_frame$boundaryNoise, 
              exp_noise=simulation_frame$experimentalNoise))
  
}

apply_perceptual_noise <- function(vector, noise){
  
  return(vector + rnorm(n=length(vector), 0, noise))
  
}

apply_boundary <- function(vector, boundary, noise){
  
  boundary.n <- rnorm(n=length(vector), boundary, noise)
  return(ifelse(vector < boundary.n, 0, 1))
  
}

apply_exp_noise <- function(vector, noise){
  
  vector <- ifelse(sample(0:1, length(vector), replace=T,
                          c(1-noise, noise)),
                   1-vector, vector)
  return(factor(vector, labels=c("A","B")))
  
}

correct_responses <- function(responses, category){
  
  correct <- ifelse(responses==category, 1, 0)
  
  if (mean(correct) < 0.5){
    require(plyr)
    responses <- revalue(responses, c("A"="B", "B"="A"))
  }
  
  return(responses)
}

point2line_distance <- function(stim_coords, gradient, y_intercept){
  ## calculates minimum distance between point and line
  ## args:
  ##   gradient: gradient of line
  ##   y_intercept: y-intercept of the line
  ##   stimulus_coordinates: a nx2 vector of stimulus coordinates
  ## output:
  ##   nx1 distances from the line, -ve above line, +ve below line
  return((gradient*stim_coords[ ,1] - stim_coords[ ,2] + y_intercept)/
           sqrt(gradient^2 + 1))
}
