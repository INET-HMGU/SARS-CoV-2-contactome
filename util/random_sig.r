# calculate the random simulation result is significant or not, compare to the original observed value

library(rethinking)

sig <- function(obs, result) {
  if (obs >= median(result)) {
    freq <- 1 - (sum(result < obs)) / length(result)
    return(round(freq, digits = 3))
  } else {
    freq <- 1 - (sum(result > obs)) / length(result)
    return(round(freq, digits = 3))
  }
}

simula <- function(time, background, query, target) {
    random <- mcreplicate(
        time,
        table(sample(background, length(query)) %in% target)[2],
        mc.cores = detectCores()
    )
    if(table(is.na(random))["FALSE"] < 10000) {
        random[is.na(random)] <- 0
    }
    return(random)
}
