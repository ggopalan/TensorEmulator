\name{tensor_emulator_comb}
\alias{tensor_emulator_comb}
\title{Emulate a simulator function using a combination of pre-trained Gaussian process regressions and random forests.}
\description{
  The function uses the pre-trained Gaussian process regressions (as output by the learn_U_GP function) and pre-trained random forests (as output by the learn_U_RF function) to emulate a simulator.}
\usage{
tensor_emulator_GP(GPs,RFs,Z,r,comb,input)
}
\arguments{
 \item{GPs}{Pre-trained GPs -- should be output of lean_U_GP.}
 \item{RFs}{Pre-trained RFs -- should be output of learn_U_RF.}
 \item{Z}{Core tensor from the truncated HOSVD. Accesible with $Z from the returned object of tensor_decomp.}
 \item{r}{This is a vector that contains the ranks used for the truncated HOSVD. Its length should be the same as the number of modes of the orignal tensor decomposed.}
 \item{comb}{This is a vector whose length is the number of modes of the tensor decomposed. If a particular element is "1", then a random forest will be used for the corresponding mode. If a particular element is "2", then a Gaussian process regression will be used for the corresponding mode.}
 \item{input}{This is a list where each item of the list corresponds to the parameter/spatio-temporal coordinate values that the user wants the emulated function value of.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo.r.
}

