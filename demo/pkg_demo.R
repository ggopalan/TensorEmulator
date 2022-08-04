#demo of TensorEmulator package

#define a test simulator that will be emulated
test_simulator <- function(x)
{
  return(sin(x[1])*cos(x[2])*tan(x[3]))
}

set.seed(1234)

#generate input values for simulator, the output of which will be stored in a tensor
n_vals <- c(20,20,20) #specify how many values of each parameter/spatio-temporal coordinate to sample
dims <- c(1,1,1) #specify the dimension of each parameter/spatio-temporal coordinate
mins <- list(c(-1), c(-1), c(-1)) #specify the minimum values of each parameter in a list
maxes <- list (c(1), c(1), c(1)) #specify the maximum values of each parameter in a list
#use the generate_params function from the package in order to use LHS to generate parameter samples
params_list <- generate_params(n_vals, dims, mins, maxes)

#build a tensor of simulator evaluations
training_tensor <- build_tensor(n_vals, test_simulator, params_list[[1]])

#use the truncated HOSVD with ranks specified
#NOTE: see Gopalan and Wikle (2021) on a procedure to determine ranks
ranks <- c(6,6,6)
decomp <- tensor_decomp(training_tensor, ranks)

#learn U matrices of the HOSVD using random forests and Gaussian process regression
RFs <- learn_U_RF(decomp, params_list[[2]], ranks)
GPs <- learn_U_GP(decomp, params_list[[2]], ranks)

#now we are ready to run the tensor emulator. Lets use the following as a test input
input <- list(c(.2), c(-.1), c(.4))
#now run the tensor emulator using...:
#all RFs:
tensor_emulator_RF(RFs, decomp$Z, ranks, input)
#all GPs:
tensor_emulator_GP(GPs, decomp$Z, ranks, input)
#RF for the first two modes and GP for the third mode
comb <- c(1,1,2)
tensor_emulator_comb(GPs,RFs,decomp$Z,ranks,comb,input)
#compare against actual simulator output
test_simulator(c(.2, -.1, .4))
