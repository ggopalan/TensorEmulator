\name{tensor_emulator_RF_ABM}
\alias{tensor_emulator_RF_ABM}
\title{Emulate an agent-based simulator function using pre-trained random forests.}
\description{
  The function uses the pre-trained random forests (as output by the learn_U_RF function) to emulate an agent-based simulator. This function uses *only* random forests. See tensor_emulator_comb for using a combination of random forests and Gaussian processes.}
\usage{
tensor_emulator_RF(RFs,Z,r,input)
}
\arguments{
 \item{RFs}{Pre-trained RFs -- should be output of learn_U_RF}
 \item{Z}{In contrast to non-ABM version, this should be the core tensor of the HOSVD decomposition multiplied by U matrices for the first, second, and third modes. This is because we are emulating the response for all spatial, agent, and temporal coordinates.}
 \item{r}{This is a vector that contains the ranks used for the truncated HOSVD for the parameter modes. Its length should be the same as the number of modes corresponding to parameters.}
 \item{input}{This is a list where each item of the list corresponds to the parameter values that the user wants the emulated function value of.}
\author{Giri Gopalan gopalan88@gmail.com}
\examples{
# See demo/pkg_demo_ABM.r.
}

