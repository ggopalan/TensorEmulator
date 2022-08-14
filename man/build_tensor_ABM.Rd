\name{build_tensor_ABM}
\alias{build_tensor_ABM}
\title{Populate a tensor with agent-based model simulator evaluations}
\description{
  This function returns a tensor object (see rTensor package) in which the elements correspond to the evaluation of an input simulator function for every combination of input parameter coordinate values. Furthermore, the first mode corresponds to spatial indices, the second mode corresponds to agent indices, the third mode corresponds to time indices, and remaining modes correspond to parameters.}
\usage{
build_tensor_ABM(n_vals,abm_dim,simulator,params)
}
\arguments{
 \item{n_vals}{Vector where each element corresponds to the number of parameter values for the parameter modes of the tensor.}
 \item{abm_dim}{Vector of length three which specifies the number of spatial indices (1,2, or 3), agent indices, and temporal indices, in that order. Assumption is that all input to the simulator is concatenated into a single vector.}
 \item{simulator}{Agent-based model simulator function. Assumption is that the function takes a single vecot input, which contains all parameters concatenated together. The assumed return of this function is a multi-dimensional array, where the first coordinate corresponds to spatial index, second corresponds to agent index, and third coordinate corresponds to the temporal index.}
 \item{params}{Matrix where each row contains a particular combination of parameter values to run the simulator on.}
 }
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo_ABM.r.
}

