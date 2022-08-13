#an example of an agent-based model simulator from Gopalan and Wikle (2022).
#load functions necessary for agent based model
Couzin_collective_v5 <- function(alpha,rho_o,s2eta,T,N,v)
{
  #Collective movement model loosely based off of Couzin et al. 2002 (J. Theo. Biol.)
  # Model includes a region of repulsion given by distances (0,alpha) and a region
  # of orientation from [alpha,rho_o) where the agents move in the average direction
  # of agents in this region; outside of this region the agents just keep
  # the same direction that they had before (plus noise given with variance s2eta)
  # ********************************************************************
  #INPUTS
  # alpha: radius of repulsion (small and fixed) [e.g., .5]
  # rho_o: radius of orientation, must be > alpha; most important parameter [e.g., 20]
  # s2eta: variance of the additive noise on location (deviation from collectiveness) [e.g., .025]
  # T: number of time steps to consider [e.g., 100]
  # N: number of agents [e.g, 10]
  # v: speed that agents move [e.g., .3]
  #
  #OUTPUTS list with:
  # x: N x (T+1) matrix of x-coordinate locations for each agent, time
  # y: N x (T+1) matrix of y-coordinate locations for each agent, time
  #
  #************************************************************************

  #initialize
  t <- 1/T  #Size of one time step
  x= matrix(0,nrow=N,ncol=(T+1))
  y= matrix(0,nrow=N,ncol=(T+1))

  #initial locations; 2 groups
  set.seed(4)

  #initial direction unit vectors
  dstrtx = runif(N,min=-1,max=1)
  dstrty = runif(N,min=-1,max=1)
  dstrt = cbind(dstrtx,dstrty)

  #initial starting locations; two groups
  x[seq(1,(N/2)),1]=runif(N/2,min=30,max=50)
  x[seq((N/2)+1,N),1]=runif(N/2,min=50,max=70)
  y[seq(1,(N/2)),1]=runif(N/2,min=30,max=50)
  y[seq((N/2)+1,N),1]=runif(N/2,min=50,max=70)

  for(t in 1:T){
    dnew = matrix(0,nrow=N,ncol=2)
    if(t == 1){
      d = dstrt
    }
    for(i in 1:N){
      #find distance
      locI <- cbind(x[i,t],y[i,t])
      fnd <- which(i != seq(1,N))
      locNI <- cbind(x[fnd,t],y[fnd,t])
      dst <- get_dist(locI,locNI)
      fnd_alpha <- which(dst < alpha)  #individuals in zone of repulsion
      fnd_rho_o <- which((alpha <= dst) & (dst < rho_o)) #zone of orientation

      r_a <-  ((matrix(1,dim(locNI)[1],1)%*%locI) - locNI)/ (t(dst) %*% matrix(1,1,2))

      if (length(fnd_alpha) != 0){
        d_a = -1 * colSums(r_a)
        dtilde = d_a }
      else if (length(fnd_rho_o) == 1) {
        d_o = d[fnd_rho_o,]
        dtilde = d_o}
      else if (length(fnd_rho_o) > 1) {
        d_o = colSums(d[fnd_rho_o,])
        dtilde = d_o}
      else {
        dtilde = d[i,]
      }
      dnew[i,] = dtilde / sqrt(sum(dtilde^2))

      x[i,t+1] = x[i,t] + dnew[i,1]*v + rnorm(1,0,sqrt(s2eta))
      y[i,t+1] = y[i,t] + dnew[i,2]*v + rnorm(1,0,sqrt(s2eta))
    }
    d <- dnew
  }

  my_list <- list(x, y)
  names(my_list) <- c("x", "y")
  return(my_list)
}

get_dist <- function(d1_locs,d2_locs){
  # Given N1 x 2 matrix of locations and N2 x 2 matrix of
  #  locations, produce a matrix of distances between all
  #  locations. Note, these two input matrices can be
  #  one and the same. The program efficiently calculates the
  #  distances using complex arithmetic; e.g., assign the x locs as the
  #  "real" and the y locs as the the "imag"

  N1 = dim(d1_locs)[1]
  N2 = dim(d2_locs)[1]

  Cd1loc <- d1_locs[,1] + sqrt(as.complex(-1))*d1_locs[,2]
  Cd2loc <- d2_locs[,1] + sqrt(as.complex(-1))*d2_locs[,2]

  tmp1 <- as.vector(Cd1loc) %*%  t(as.vector(rep(1,N2)))
  tmp2 <- t(as.vector(Cd2loc) %*% t(as.vector(rep(1,N1))))

  DST <- Mod(tmp1 - tmp2)  #% Distance matrix
  return(DST)
}
#wrapper function for ABM simulator
ABM_simulator <- function(theta)
{
  rho <- theta[1]
  v <- theta[2]
  out <- Couzin_collective_v5(.5,rho,.025,100,20,v)
  out <- array(as.vector(rbind(as.vector(out$x), as.vector(out$y))),dim=c(2,20,101))
  return(out)
}
