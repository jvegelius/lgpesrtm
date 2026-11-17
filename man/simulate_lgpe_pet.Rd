\name{simulate_lgpe_pet}
\alias{simulate_lgpe_pet}
\title{Simulate dynamic PET data}
\description{Simulate dynamic PET data.}
\usage{
simulate_lgpe_pet(n = 50, t_obs = seq(0.5, 60.5, by = 1), onset_time = 40,
  activation_amplitude = 0.012, activation_shape = "gamma",
  common_activation = TRUE, seed = NULL)
}
\arguments{
  \item{n}{Number of subjects.}
  \item{t_obs}{Observation times in minutes.}
  \item{onset_time}{Intervention onset in minutes.}
  \item{activation_amplitude}{Amplitude of the group activation curve.}
  \item{activation_shape}{Shape of the activation curve, such as gamma or logistic.}
  \item{common_activation}{Logical flag controlling whether activation is shared within the activation group.}
  \item{seed}{Optional random seed.}
}
\value{A list with simulated data, truth and the inputs for fitting.}
