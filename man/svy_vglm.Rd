\name{svy_vglm}
\alias{svy_vglm}
\alias{svy_vglm.survey.design}
\alias{svy_vglm.svyrep.design}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design-based inference for vector generalised linear models
}
\description{
This function provides design-based (survey) inference for Thomas Yee's
vector generalised linear models.  It works by calling \code{vglm} with
sampling weights, and then either using resampling (replicate weights)
or extracting the influence functions and using a Horvitz-Thompson-type
sandwich estimator.

}
\usage{
svy_vglm(formula, family, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula, as for \code{\link{vglm}}
  }
  \item{family}{
Model family, as for \code{\link{vglm}}
  }
  \item{design}{
Survey design object
  }
  \item{\dots}{
Other arguments to pass to  \code{\link{vglm}}
  }
}

\value{
An S3 object of class \code{svy_glm} with \code{print}, \code{coef} and \code{vcov}
methods, containing the design in the \code{design} component and a
fitted \code{vglm} object in the \code{fit} component.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Thomas Lumley
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, family, design, ...) 
{
    UseMethod("svy_vglm", design)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
