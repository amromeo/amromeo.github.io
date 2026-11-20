\name{simulate_ri_study}
\alias{simulate_ri_study}
\alias{simulate_ri_case}
\title{Simulate paired clinic and ward studies with a known reference interval}
\usage{
simulate_ri_study(scenario = 3L, replicate = 1L, n_clinic = 250L,
                  n_ward = 750L, n_anchor = 40L,
                  selection = c("clinical", "in_range"))
simulate_ri_case(scenario, replicate, n_clinic = 1000L, n_ward = 3000L)
}
\arguments{
  \item{scenario}{1 through 4. See Details.}
  \item{replicate}{Positive integer draw identifier used to derive the seeds.}
  \item{n_clinic, n_ward}{Numbers of independent routine results per setting.}
  \item{n_anchor}{Number of additional reference individuals, including zero.}
  \item{selection}{\code{"clinical"} draws screened values from the reference population. \code{"in_range"} also restricts reported values to 26 through 34 seconds, simulating selection bias.}
}
\description{Creates synthetic aPTT-like results in seconds. Reference values follow a normal distribution with mean 30 and SD 3. Results are rounded to one decimal place. The target is the central 95 percent interval of the underlying continuous reference distribution. The generator preserves the caller's random number state.}
\details{Scenario 1 has affected fractions of 5 and 10 percent in clinic and ward, with affected mean 42 and SD 4. Scenario 2 has 20 and 75 percent affected, with affected mean 36 and SD 4.5. Scenario 3 has the same fractions, with 55 percent of affected ward results drawn from a lower component (mean 24, SD 2.5). Scenario 4 has 40 percent affected in both settings and the same upper component as scenario 2. Screened draws are independent of routine results. At a fixed scenario and replicate, smaller screened groups use a prefix of the same draw.}
\value{\code{simulate_ri_study} returns \code{routine} (a data frame with \code{value} and \code{setting}), \code{anchor} and \code{truth}. \code{simulate_ri_case} returns the routine data frame alone.}
\examples{
d <- simulate_ri_study()
table(d$routine$setting)
d$truth
}
