\name{spT.priors}
\alias{spT.priors}
\title{Priors for the spatio-temporal models.}
\description{This command is useful to assign the hyper-parameters of the prior distributions.
}
\usage{
spT.priors(model="GP", inv.var.prior=Gamm(a=2,b=1),
  beta.prior=Norm(0,10^10), rho.prior=Norm(0,10^10))
}
\arguments{
\item{model}{The spatio-temporal models, current input: "GP", "AR", and "GPP".}
\item{inv.var.prior}{The hyper-parameter for the Gamma prior distribution (with mean = a/b) of the precision (inverse variance) model parameters (e.g., 1/\eqn{\sigma}2_\eqn{\epsilon}, 1/\eqn{\sigma}2_\eqn{\eta}).}
\item{beta.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\beta} model parameters.}
\item{rho.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\rho} model parameter.}
}
\note{
If no prior information are given (assigned as NULL), then it use flat prior values of the corresponding distributions. \cr
\code{Gamm} and \code{Norm} refers to Gamma and Normal distributions respectively.
}
\seealso{
\code{\link{spT.Gibbs}, \link{predict.spT}, \link{spT.initials}}.
}
\examples{
\dontrun{
## 

priors<-spT.priors(model="GPP",inv.var.prior=Gamm(2,1),
      beta.prior=Norm(0,10^4))
priors

##
}
}
\keyword{spT}
