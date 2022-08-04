\name{tensor_emulator_RF}
\alias{tensor_emulator_RF}
\title{Emulate a simulator function using pre-trained random forests.}
\description{
  The function uses the pre-trained random forests (as output by the learn_U_RF function) to emulate a simulator. This function uses *only* random forests. See tensor_emulator_comb for using a combination of random forests and Gaussian processes.}
\usage{
tensor_emulator_RF(RFs,Z,r,input)
}
\arguments{
 \item{RFs}{Pre-trained RFs -- should be output of learn_U_RF}
 \item{Z}{Core tensor from the truncated HOSVD. Accesible with $Z from the returned object of tensor_decomp.}
 \item{r}{This is a vector that contains the ranks used for the truncated HOSVD. Its length should be the same as the number of modes of the orignal tensor decomposed.}
 \item{input}{This is a list where each item of the list corresponds to the parameter/spatio-temporal coordinate values that the user wants the emulated function value of.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

