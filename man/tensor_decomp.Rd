\name{tensor_decomp}
\alias{tensor_decomp}
\title{Perform a truncated higher-order singular value decomposition}
\description{
  This function returns the truncated higher-order singular value decomposition on the input tensor with ranks specified by user. This is a  wrapper function for hosvd from rTensor.}
\usage{
tensor_decomp(training_tensor, r)
}
\arguments{
 \item{training_tensor}{Tensor object (see rTensor) that stores simulator evaluations, ideally the output of build_tensor.}
 \item{r}{Vector that specifies the ranks to be used in the truncated higher-order singular value decomposition. Length of vector should match with the number of modes of training_tensor.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

