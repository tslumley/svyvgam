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

\seealso{
\code{\link{nhanes_sxq}}

\code{\link{vglm}}

\code{\link{svydesign}}
\code{\link{svrepdesign}}


}

\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

## Ordinary Gaussian regression
m1<-svyglm(api00~api99+mobility+ell, design=dclus2,family=gaussian)
## same model, but with the variance as a second parameter
m2<-svy_vglm(api00~api99+mobility+ell, design=dclus2,family=uninormal())
m1
m2
SE(m1)
SE(m2)

summary(m1)
summary(m2)

## Proportional odds model
dclus2<-update(dclus2, mealcat=as.ordered(cut(meals,c(0,25,50,75,100))))
a<-svyolr(mealcat~avg.ed+mobility+stype, design=dclus2)
b<-svy_vglm(mealcat~avg.ed+mobility+stype, design=dclus2, family=propodds())
a
b
SE(a)
SE(b) #not identical, because svyolr() uses approximate Hessian


## Zero-inflated Poisson
data(nhanes_sxq)
nhdes = svydesign(id=~SDMVPSU,strat=~SDMVSTRA,weights=~WTINT2YR,
       nest=TRUE, data=nhanes_sxq)

sv1<-svy_vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC,
      zipoisson(), design=nhdes, crit = "coef")
sv1
summary(sv1)

## Multinomial
## Reference group (non-Hispanic White) average older and more educated
## so coefficients are negative
mult_eth<- svy_vglm(RIDRETH1~RIDAGEYR+DMDEDUC,
     family=multinomial(refLevel=3), design=nhdes)

## separate logistic regressions are close but not identical
two_eth<-svyglm(I(RIDRETH1==1)~RIDAGEYR+DMDEDUC, family=quasibinomial,
    design=subset(nhdes, RIDRETH1 \%in\% c(1,3)))

summary(mult_eth)
summary(two_eth)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
