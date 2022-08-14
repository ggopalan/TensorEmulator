#TensorEmulator package functions

#build a training tensor
#NOTE: simulator should take a single vector as input. The ordering of this vector should
#correspond to the ordering of elements in params.
#For example: (x,y,z,t,theta_11,theta_12,theta_2).
build_tensor <- function(n_vals, simulator, params)
{
  print("Evaluating simulator...")
  simulator_evals <- apply(params,1,simulator)
  training_tensor <- array(simulator_evals,dim=n_vals)
  return(as.tensor(training_tensor))
}

#function to apply tensor decomposition using hosvd with ranks R
tensor_decomp <- function(training_tensor,r)
{
  return(hosvd(training_tensor,ranks = r))
}

#learn U functions with random forests
learn_U_RF <- function(decomp, params_list, r)
{
  K <- length(params_list)
  RForests <- vector(mode="list", length= K)
  for(i in 1:K)
  {
    RForests[[i]] <- vector(mode = "list", length = r[i])
    for(index in 1:r[i])
    {
      RForests[[i]][[index]] <- randomForest(params_list[[i]],decomp$U[[i]][,index],ntree=1000)
    }
  }
  return(RForests)
}

#learn U functions with GPs
learn_U_GP <- function(decomp, params_list, r)
{
  K <- length(params_list)
  GPs <- vector(mode="list", length= K)
  for(i in 1:K)
  {
    GPs[[i]] <- vector(mode = "list", length = r[i])
    for(index in 1:r[i])
    {
      GPs[[i]][[index]] <- gausspr(x=params_list[[i]],y = decomp$U[[i]][,index])
    }
  }
  return(GPs)
}

#emulator function using pure RFs
tensor_emulator_RF <- function(RFs,Z,r,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
    cur_val <- input[[i]]
    cur_u <- rep(0,r[[i]])
    for(j in 1:r[[i]])
    {
      cur_RF <- RFs[[i]][[j]]
      cur_u[j] <- predict(cur_RF,cur_val)
    }
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i)
  }
  return(predict_tensor)
}

#emulator function using pure GPs
tensor_emulator_GP <- function(GPs,Z,r,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
    cur_val <- input[[i]]
    cur_u <- rep(0,r[[i]])
    for(j in 1:r[[i]])
    {
      cur_GP <- GPs[[i]][[j]]
      cur_u[j] <- predict(cur_GP,cur_val)
    }
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i)
  }
  return(predict_tensor)
}

#emulator function using a user-specified combination of RFs and GPs
tensor_emulator_comb <- function(GPs,RFs,Z,r,comb,input)
{
  predict_tensor <- Z
  K <- length(r)
  for(i in 1:K)
  {
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
    predict_tensor <- ttm(predict_tensor, mat=matrix(cur_u,nrow=1),m=i)
  }
  return(predict_tensor)
}


#Function that will generate spatio-temporal values and parameters
#using Latin hypercube sampling for each mode, to populate the tensor
generate_params <- function(n_vals,dims,mins,maxes)
{
  K <- length(n_vals)
  params_list <-vector(mode="list",length=length(n_vals))
  print("Performing LHS for each tensor mode")
  for(i in 1:K)
  {
    cur_min <- mins[[i]]
    cur_max <- maxes[[i]]
    total_length <- cur_max-cur_min
    #use LHS to generate values
    cur_vals <- randomLHS(n_vals[i],dims[i])
    for(j in 1:length(total_length))
      {
        cur_vals[,j] <- total_length[j]*cur_vals[,j]
        cur_vals[,j] <- cur_vals[,j]+cur_min[j]
      }
    params_list[[i]] <- cur_vals
  }
  nval_list <-vector(mode="list",length=length(n_vals))
  for(i in 1:K)
  {
    nval_list[[i]] <- seq(n_vals[i])
  }
  #now create one large matrix where each row is a different combination
  #relies on expand.grid
  print("Populating param matrix")
  params <- matrix(rep(0,sum(dims)*prod(n_vals)), ncol = sum(dims))
  index <- 1
  indices <- expand.grid(nval_list)
  for(index in 1:(dim(indices)[1]))
  {
    cur_comb <- c()
    for(i in 1:K)
    {
      cur_comb <- c(cur_comb,params_list[[i]][indices[index,i],])
    }
    params[index,] <- cur_comb
  }
  return(list(params, params_list))
}


