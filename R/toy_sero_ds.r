##' Create a toy individualised data set of seroprevalence
##'
##' @param n number of sample to generate per country
##' @param pos_eqi a named list of vectors of size 3, each name corresponding to the country represented, and the three elements of the vector corresponding to the proportion positive and equivocal (the reminder being negative)
##' @param survey.yr an integer vector: year(s) during which the seroprevalence studies occurred
##' @param lower.age.limits an integer vector: lower limits of the age groups for which to generate data
##' @return a data frame representing a toy data set of individual serology results, with columns \code{survey.yr} (the year of the survey), \code{country} (the country), \code{age1} (the lower limit of the age group) and \code{stdres} (the standardised result)
##' @author Sebastian Funk <sebastian.funk@lshtm.ac.uk>
##' @export
##' @examples
##'   ms.sero <- toy_sero_ds(survey.yr=2002, pos_eqi=list(Spain=c(0.8, 0.1), Ireland=c(0.9, 0)), lower.age.limits=seq(0, 70, by=5), n=100)
toy_sero_ds <- function(n, pos_eqi, survey.yr, lower.age.limits)
{
  if (length(survey.yr)==1) survey.yr <- rep(survey.yr, length(pos_eqi))

  df <- list()
  for (i in seq_len(length((pos_eqi))))
  {
      country <- names(pos_eqi)[i]
      df[[i]] <-
          data.frame(survey.yr = survey.yr[i],
                     country = country,
                     age1 = sample(lower.age.limits,
                                   size=n,
                                   replace=TRUE), 
                     stdres = sample(x=c("POS", "EQI", "NEG"),
                                     size=n,
                                     prob=c(pos_eqi[[country]], 1-sum(pos_eqi[[country]])),
                                     replace=TRUE))
  }

  return(do.call(rbind, df))
}
