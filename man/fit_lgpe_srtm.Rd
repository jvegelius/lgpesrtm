\name{fit_lgpe_srtm}
\alias{fit_lgpe_srtm}
\title{Fit the LGPE-SRTM model}
\description{Fit the LGPE-SRTM model.}
\usage{
fit_lgpe_srtm(df, predvars, par_in, estpat, control)
}
\arguments{
  \item{df}{Data frame of observed PET curves with subject identifiers, time, onset, and measurement information.}
  \item{predvars}{List of predictor names for R1, k2 and k2a, respectively.}
  \item{par_in}{List of initial model parameters.}
  \item{estpat}{List of estimation restrictions.}
  \item{control}{List of fitting controls.}
}
\value{A list of fitted kinetic parameters, uncertainty estimates and test results.}
