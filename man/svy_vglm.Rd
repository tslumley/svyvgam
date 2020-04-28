\name{svy_vglm}
\alias{svy_vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design-based inference for vector generalised linear models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
