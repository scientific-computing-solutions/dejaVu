#' Simulated recurrent event data.
#'
#' A simulated dataset containing a randomised treatment group, follow-up time,
#' and number of events, for 500 patients. The planned follow-up period for the 
#' study was 1 year, but some patients dropped out early and so their
#' follow-up ended prematurely (i.e. before 1 year)
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{z}{a binary variable indicating randomised treatment group}
#'   \item{y}{number of events observed during patient's follow-up}
#'   \item{fupTime}{the time in years the patient was followed up for}
#'   ...
#' }
#' @source Simulated data
"simData"
