\name{generate_params}
\alias{generate_params}
\title{Sample parameters/spatio-temporal coordinates}
\description{
  This function samples parameters/spatio-temporal coordinates using Latin hypercube sampling with the lhs package. It returns both a matrix that has every combination of parameter/spatio-temporal coordinates in each row, and a list of matrices for each parameter/spatio-temporal coordinate sampled. Both the matrix and list are further used by other functions of the package.}
\usage{
generate_params(n_vals, dims, mins, maxes)
}
\arguments{
 \item{n_vals}{Vector where each element corresponds to the number of parameter/spatio-temporal coordinate values for each mode of the tensor. That is, specify how many parameters/spatio-temporal coordinates to sample using latin hyper-cube sampling.}
 \item{dims}{Vector where each elements corresponds to the dimension of the corresponding parameter/spatio-temporal coordinate. For example, if there are three modes with the first mode a vector of length 2, and the second and third modes scalars, then dim should be c(2,1,1).}
 \item{mins}{List where each component corresponds to the minimum parameter/spatio-temporal coordinate value. Needed to specify the range of sampling.}
 \item{maxes}{List where each corresponds to the maximum parameter/spatio-temporal coordinate value. Needed to the specify the range of sampling.}
 }
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

