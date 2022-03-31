# Required package

library(plyr)

##################################################

# Variable declaration

theta0 <- 0
theta1 <- angle(c(1, 0), c(2, 1))
theta2 <- angle(c(2, 1), c(1, 1))
theta3 <- angle(c(1, 1), c(-1, 0))
theta4 <- angle(c(-1, 0), c(-2, -1))
theta5 <- angle(c(-2, -1), c(-1, -1))
theta6 <- angle(c(-1, -1), c(1, 0))
ang0 <- 0
ang1 <- theta1
ang2 <- ang1 + theta2
ang3 <- ang2 + theta3
ang4 <- ang3 + theta4
ang5 <- ang4 + theta5
ang6 <- ang5 + theta6

##################################################

# Function definition


GPT <- function(signal){
  
  n <- length(signal) - 2 # amount of points after applying the GPT
  
  points <- matrix(0, nrow=n, ncol=2) # initialization
  
  for (i in 1:n){
    points[i,1] <- 2*signal[i+2]-signal[i+1]-signal[i] # first component in the plane
    points[i,2] <- signal[i+2]-signal[i] # second component in the plane
  }
  
  return(points) 
}


cartesian2polar <- function(x, y){
  
  # Radius
  r <- sqrt(x^2 + y^2)
  
  # Origin
  if (x == 0 & y == 0){
    t <- 0
  }
  
  # Angle of a point over the posive y axis
  if (x == 0 & y > 0){
    t <- pi / 2
  }
  
  # Angle of a point over the negative y axis
  if (x == 0 & y < 0){
    t <- 3 * pi / 2
  }
  
  # Angle of a point over the posive x axis
  if (y == 0 & x > 0){
    t <- 0
  }
  
  # Angle of a point over the negative x axis
  if (y == 0 & x < 0){
    t <- pi
  }
  
  # Angle of a point in the first quadrant
  if (x > 0 & y > 0){
    t <- atan(y / x)
  }
  
  # Angle of a point in the second quadrant
  if (x < 0 & y > 0){
    t <- pi + atan(y / x)
  }
  
  # Point in the third quadrant
  if (x < 0 & y < 0){
    t <- pi + atan(y / x)
  }
  
  # Angle of a point in the fourth quadrant
  if (x > 0 & y < 0){
    t <- 2*pi + atan(y / x)
  }
  
  return(c(r,t))
}


vec_norm <- function(v){
  
  return(sqrt(v[1] * v[1] + v[2] * v[2])) 
  
}


angle <- function(v, w){
  
  cdot <- (v[1] * w[1] + v[2] * w[2]) # inner product
  vnorm <- vec_norm(v) # norm of v
  wnorm <- vec_norm(w) # norm of w
  
  theta <- acos(cdot / (vnorm * wnorm)) # angle bewteen v and w
  
  return(theta)
}



pattern <- function(ang){
  
  if (ang0 <= ang & ang < ang1){
    pat <- 213
  }
  
  if (ang1 <= ang & ang < ang2){
    pat <- 123
  }
  
  if (ang2 <= ang & ang < ang3){
    pat <- 132
  }
  
  if (ang3 <= ang & ang < ang4){
    pat <- 231
  }
  
  if (ang4 <= ang & ang < ang5){
    pat <- 321
  }
  
  if (ang5 <= ang & ang < ang6){
    pat <- 312
  }
  
  return(pat)
}


partition <- function(set, n_part){
  
  m <- min(set) # minimum of the set
  M <- max(set) # maximum of the set
  part_length <- (M-m)/n_part # part length
  
  separation <- seq(m, M, part_length) # partition limits
  
  q <- length(separation) # number of parts
  
  part <- vector() # initialization
  
  # Partition assignment
  for (i in 1:length(set)){
    
    if (set[i] == m){
      part[i] <- 1
    } else if (set[i] == M){
      part[i] <- n_part
    } else {
      sorted <- sort(c(set[i], separation))
      if (length(sorted) == q + 1){
        part[i] <- min(which(sorted == set[i])) - 1
      } else {
        part[i] <- which(sep == set[i]) - 1
      }
    }
  }
  
  return(part)
}


partition_length <- function(set){
  
  nset <- length(set) # number of elements in the set
  opt_length <- 1 # initialization
  
  repeat{
    pp <- partition(set, opt_length)
    if (length(unique(pp)) == opt_length){
      opt_length <- opt_length + 1
    } else {
      break
    } 
  }
  
  return(opt_length-1)
  
}


voe <- function(signal){
  
  # GPT points
  gpt <- GPT(signal)
  
  # Patterns assignment
  gpt_pat <- vector()
  for (i in 1:dim(gpt)[1]){
    polar_coord <- cartesian2polar(gpt[i, 1], gpt[i, 2])
    gpt_pat[i] <- pattern(polar_coord[2])
  }
  
  # Data arrangement
  data_points <- as.data.frame(cbind(gpt, gpt_pat))
  colnames(data_points) <- c("Delta_M", "Delta_S", "Pattern")
  
  # Pattern splittings
  by_pattern <- split(data_points, data_points$Pattern)
  m <- length(by_pattern)
  
  # Vector of entropy computation
  vec_entropy <- rep(NA, 6)
  
  for (i in 1:m){
    
    # Data preparation
    data_points_splittings <- as.data.frame(by_pattern[i])
    colnames(data_points_splittings) <- colnames(data_points)
    
    if (data_points_splittings$Pattern[1] == 213){
      vec_coord <- 1
    } else if (data_points_splittings$Pattern[1] == 123){
      vec_coord <- 2
    } else if (data_points_splittings$Pattern[1] == 132){
      vec_coord <- 3
    } else if (data_points_splittings$Pattern[1] == 231){
      vec_coord <- 4
    } else if (data_points_splittings$Pattern[1] == 321){
      vec_coord <- 5
    } else if (data_points_splittings$Pattern[1] == 312){
      vec_coord <- 6
    }
    
    # Partition of the set of norms
    number_points <- dim(data_points_splittings)[1]
    
    norms <- vector()
    for (j in 1:number_points){
      norms[j] <- vec_norm(data_points_splittings[j, 1:2])
    }
    norms <- unlist(norms)
    
    number_part <- partition_length(norms)
    parts <- partition(norms, number_part)
    
    data_points_splittings <- cbind(data_points_splittings, parts)
    
    # Probability function
    frequencies <- count(data_points_splittings, vars = "parts")["freq"]
    prob <- frequencies / sum(frequencies)
    
    # Entropy
    H <- -sum(prob * log(prob))
    if (H > 0){
      vec_entropy[vec_coord] <- H / log(dim(prob)[1])
    } else {
      vec_entropy[vec_coord] <- 0
    }
  }
  
  return(vec_entropy)
}