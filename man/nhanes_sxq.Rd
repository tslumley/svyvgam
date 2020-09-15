\name{nhanes_sxq}
\alias{nhanes_sxq}
\docType{data}
\title{
Data from NHANES: number of sex partners
}
\description{
These data are from the NHANES 2003-2004 survey in the US. They provide an example of  overdispersed count data that motivates a two-component zero-inflation model

}
\usage{data("nhanes_sxq")}
\format{
  A data frame with 2992 observations on the following 7 variables.
  \describe{
    \item{\code{SDMVPSU}}{Primary Sampling Unit}
    \item{\code{SDMVSTRA}}{stratum}
    \item{\code{WTINT2YR}}{weights}
    \item{\code{malepartners}}{lifetime number of male sexual partners}
    \item{\code{RIDAGEYR}}{age in years}
    \item{\code{DMDEDUC}}{level of  education: 1=less than high school, 2=high school, 3-more than high school, 7=refused}
    \item{\code{RIDRETH1}}{Race/ethnicity: 1=Mexican American, 2=Other Hispanic, 4=non-Hispanic White, 5=non-Hispanic Black, 5=Other}
  }
}

\source{
NHANES files \code{demo_c.xpt} and \code{sxq_c.xpt}
}
\seealso{
Construction of the data set is described by \url{https://notstatschat.rbind.io/2015/05/26/zero-inflated-poisson-from-complex-samples/}

}

\examples{
data(nhanes_sxq)
nhdes = svydesign(id=~SDMVPSU,strat=~SDMVSTRA,weights=~WTINT2YR, nest=TRUE, data=nhanes_sxq)
svy_vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), design=nhdes, crit = "coef")
}
\keyword{datasets}
