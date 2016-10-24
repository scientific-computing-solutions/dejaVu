##' Expand event counts into a list of event times
##'
##' This function exists to allow clinical trial data which typically
##' gives event counts over time to be plugged into this software, which
##' relies on actual event counts.
##'
##' This function always produces a warning: anyone relying on this
##' function to actually analyze data should take great care.
##'
##' @param count a vector of event counts.  All entries must be
##' non-negative.
##' @param time a matching (strictly positive) vector of followup
##' times.
##' @examples
##' expandEventCount(count=c(0, 20), time=c(10, 20))
##' @return a list of vectors of event times
##' @export
expandEventCount <- function(count, time) {
    if ( (!is.numeric(count)) || (!is.numeric(time)) ) {
        stop("Require numeric input!")
    }
    true.count <- as.integer(count)
    ## conversion to integer is exact up to any number we're going to
    ## use here
    if (any(as.numeric(true.count) != count)) {
        stop("Event counts are non-integer.")
    }
    if (any(true.count < 0)) {
        stop("Cannot have negative event counts")
    }
 
    time <- as.numeric(time)
    if (any(time <= 0)) {
        stop("Cannot have non-positive followup times")
    }
    if (length(time) == 1) {
        time <- rep.int(time, length(true.count))
    }
    if (length(time) != length(true.count)) {
        stop("time does not match event count")
    }

    warning("Synthesizing fake event times from event count data.")
    
    mapply(function(count, duration) {
               duration * (seq_len(count) / count)
           },
           true.count, time, SIMPLIFY=FALSE)
}

