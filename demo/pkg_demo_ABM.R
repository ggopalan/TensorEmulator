#demo of TensorEmulator package
#agent based model emulation

# use R package functions to build training tensor
n_vals <- c(2^3,2^3) #how many parameter values should be sampled
dims <- c(1,1) #specify the dimension of each parameter
mins <- list(c(1), c(.1)) #specify the minimum values of each parameter in a list
maxes <- list (c(40),(1)) #specify the maximum values of each parameter in a list
params_list <- generate_params(n_vals, dims, mins, maxes)
abm_dim <- c(2,20,101) #specify dimensions of the ABM output (space,agents,time)
#build training tensor using ABM version
training_tensor <- build_tensor_ABM(n_vals,abm_dim,ABM_simulator, params_list[[1]])
ranks <- c(5,5)
full_ranks <- c(abm_dim,ranks)
#use HOSVD to decompose tensor
full_decomp <- tensor_decomp(training_tensor, full_ranks)
#extract U matrices for parameter components of decomposition
U <-  vector(mode = "list", length = length(ranks))
for(i in 1:length(ranks))
{
  U[[i]] <- full_decomp$U[[length(abm_dim)+i]]
}
decomp <- list(U)
names(decomp) <- c('U')
#learn U with RF
#learn U matrices of the HOSVD using random forests and Gaussian process regression
RFs <- learn_U_RF(decomp, params_list[[2]], ranks)
GPs <- learn_U_GP(decomp, params_list[[2]], ranks)

#now we are ready to run the tensor emulator. Lets use the following as a test input
input <- list(c(35), c(.5))
Z <- full_decomp$Z
Z <- ttm(Z, mat = full_decomp$U[[1]], m = 1)
Z <- ttm(Z, mat = full_decomp$U[[2]], m = 2)
Z <- ttm(Z, mat = full_decomp$U[[3]], m = 3)
#now run the tensor emulator using...:
#all RFs:
tensor_emulator_RF_ABM(RFs, Z, ranks, input)
#all GPs:
tensor_emulator_GP_ABM(GPs, Z, ranks, input)
#RF for the first  modes and GP for the second mode
comb <- c(1,2)
tensor_emulator_comb_ABM(GPs,RFs,Z,ranks,comb,input)
#compare against actual simulator output
ABM_simulator(c(35,.5))
