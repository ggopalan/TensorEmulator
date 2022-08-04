\name{tensor_emulator_GP}
\alias{tensor_emulator_GP}
\title{Emulate a simulator function using pre-trained Gaussian process regressions.}
\description{
  The function uses the pre-trained Gaussian process regressions (as output by the learn_U_GP function) to emulate a simulator. This function uses *only* Gaussian process regressions. See tensor_emulator_comb for using a combination of random forests and Gaussian processes.}
\usage{
tensor_emulator_GP(GPs,Z,r,input)
}
\arguments{
 \item{GPs}{Pre-trained GPs -- should be output of learn_U_GP.}
 \item{Z}{Core tensor from the truncated HOSVD. Accesible with $Z from the returned object of tensor_decomp.}
 \item{r}{This is a vector that contains the ranks used for the truncated HOSVD. Its length should be the same as the number of modes of the orignal tensor decomposed.}
 \item{input}{This is a list where each item of the list corresponds to the parameter/spatio-temporal coordinate values that the user wants the emulated function value of.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

