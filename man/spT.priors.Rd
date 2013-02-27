\name{spT.priors}
\alias{spT.priors}
\title{Priors for the spatio-temporal models.}
\description{This command is useful to assign the hyper-parameters of the prior distributions.
}
\usage{
spT.priors(model,var.prior=Gam(a=2,b=1),
  beta.prior=Nor(0,10^10), rho.prior=Nor(0,10^10),
  phi.prior=Gam(a=2,b=1))
}
\arguments{
\item{model}{The spatio-temporal models, current input: "GP", "AR", and "GPP".}
\item{var.prior}{The hyper-parameter for the Gamma prior distribution (with mean = a/b) of the variance model parameters (e.g., \eqn{\sigma}2_\eqn{\epsilon}, \eqn{\sigma}2_\eqn{\eta}).}
\item{beta.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\beta} model parameters.}
\item{rho.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\rho} model parameter.}
\item{phi.prior}{The hyper-parameter for the Gamma prior distribution (with mean = a/b) of the spatial decay (\eqn{\phi}) parameter.}
}
\note{
If no prior information are given (assigned as NULL), then it use flat prior values of the corresponding distributions. \cr
\code{Gam} and \code{Nor} refers to Gamma and Normal distributions respectively.
}
\seealso{
\code{\link{spT.Gibbs}, \link{predict.spT}, \link{spT.prediction}, \link{spT.forecast}, \link{spT.initials}}.
}
\examples{
\dontrun{
## 

priors<-spT.priors(model="GPP",var.prior=Gam(2,1),
      beta.prior=Nor(0,10^4))
priors

##
}
}
\keyword{spT}
