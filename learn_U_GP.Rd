\name{learn_U_GP}
\alias{learn_U_GP}
\title{Use Gaussian process regression to learn U matrices of the truncated higher-order singular value decomposition}
\description{
  This function uses the kernlab R package to learn the U matrices from the truncated HOSVD, with the gausspr command. The returned object is a list of lists -- the length of the first list is the number of modes of the training_tensor (number of U matrices), whereas the second list has a length based on the rank of the U matrix.}
\usage{
learn_U_GP(decomp, params_list,r)
}
\arguments{
 \item{training_tensor}{Tensor truncated HOSVD decomposition as returned by tensor_decomp.}
 \item{params_list}{This is a list whose length is the number of modes of the original tensor that has been decomposed. Each element of the list is a matrix which has each parameter/spatio-temporal coordinate sampled on a seperate row.}
 \item{r}{This is a vector that contains the ranks used for the truncated HOSVD. Its length should be the same as the number of modes of the orignal tensor decomposed.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}
