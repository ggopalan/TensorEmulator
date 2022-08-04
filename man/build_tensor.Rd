\name{build_tensor}
\alias{build_tensor}
\title{Populate a tensor with simulator evaluations}
\description{
  This function returns a tensor object (see rTensor package) in which the elements correspond to the evaluation of an input simulator function for every combination of input parameter/spatio-temporal coordinate values.}
\usage{
build_tensor(n_vals, simulator, params)
}
\arguments{
 \item{n_vals}{Vector where each element corresponds to the number of parameter/spatio-temporal coordinate values for each mode of the tensor.}
 \item{simulator}{Simulator function that is to be emulated. Assumption is that all input to the simulator is concatenated into a single vector.}
 \item{params}{Matrix where each row contains a particular combination of parameter/spatio-temporal values to run the simulator on.}
 }
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

