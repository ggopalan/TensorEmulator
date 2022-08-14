#functions specific to emulating agent-based model (ABM) simulators

#Another function to populate a tensor specifically with agent-based model simulations.
#Assumption is that simulator function returns a multidimensional array and takes a single
#vector as input. First dimension of returned array is number of spatial coordinates (1, 2, or 3),
#second dimension is number of agents, third dimension is number of time indices.
build_tensor_ABM <- function(n_vals,abm_dim,simulator,params)
{
  print("Evaluating simulator...")
  simulator_evals <- as.vector(apply(params,1,simulator))
  training_tensor <- array(simulator_evals,dim=c(abm_dim,n_vals))
  return(as.tensor(training_tensor))
}

#emulator function using pure RFs
tensor_emulator_RF_ABM <- function(RFs,Z,r,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
    print(i)
    cur_val <- input[[i]]
    cur_u <- rep(0,r[[i]])
    for(j in 1:r[[i]])
    {
      cur_RF <- RFs[[i]][[j]]
      cur_u[j] <- predict(cur_RF,cur_val)
    }
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i+3)
  }
  return(predict_tensor)
}

#emulator function using pure GPs
tensor_emulator_GP_ABM <- function(GPs,Z,r,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
    print(i)
    cur_val <- input[[i]]
    cur_u <- rep(0,r[[i]])
    for(j in 1:r[[i]])
    {
      cur_GP <- GPs[[i]][[j]]
      cur_u[j] <- predict(cur_GP,cur_val)
    }
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i+3)
  }
  return(predict_tensor)
}

#emulator function using a user-specified combination of RFs and GPs
tensor_emulator_comb_ABM <- function(GPs,RFs,Z,r,comb,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
    print(i)
    cur_val <- input[[i]]
    cur_u <- rep(0,r[[i]])
    #specify RF or GP via comb
    if(comb[i] == 1)
    {
      cur_model_type <- RFs[[i]]
    }
    else
    {
      cur_model_type <- GPs[[i]]
    }
    for(j in 1:r[[i]])
    {
      cur_model <- cur_model_type[[j]]
      cur_u[j] <- predict(cur_model,cur_val)
    }
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i+3)
  }
  return(predict_tensor)
}
