\name{fit_lpntpet_catalog}
\alias{fit_lpntpet_catalog}
\title{Fit lp-ntPET using a basis catalog}
\description{Fit lp-ntPET using a basis catalog.}
\usage{
fit_lpntpet_catalog(time, ct, cr, t_onset, weights = NULL, basis_catalog = NULL)
}
\arguments{
  \item{time}{Observation times.}
  \item{ct}{Target time-activity curve.}
  \item{cr}{Reference time-activity curve.}
  \item{t_onset}{Intervention onset time.}
  \item{weights}{Optional regression weights.}
  \item{basis_catalog}{Optional precomputed catalog.}
}
\value{A list with estimates from the selected basis fit.}
