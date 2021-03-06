##' Simulate chains using a branching process
##'
##' @param n number of simulations to run.
##' @param offspring offspring distribution: a character string corresponding to
##'   the R distribution function (e.g., "pois" for Poisson, where
##'   \code{\link{rpois}} is the R function to generate Poisson random numbers) 
##' @param stat statistic to calculate ("size" or "length" of chains)
##' @param infinite a size or length from which the size/length is to be
##'     considered infinite
##' @param tree return the tree of infectors
##' @param serial the serial interval; a function that takes one parameter
##' (`n`), the number of serial intervals to randomly sample; if this parameter
##'   is set, `chain_sim` returns times of infection, too; implies (`tree`=TRUE)
##' @param t0 start time (if serial interval is given); either a single value (0
##'     by default for all simulations, or a vector of length `n` with initial
##'     times) 
##' @param tf end time (if serial interval is given)
##' @param ... parameters of the offspring distribution
##' @return a vector of sizes/lengths (if \code{tree==FALSE} and no serial
##'   interval given), or a data frame with columns `n` (simulation ID), `time`
##'   (if the serial interval is given) and (if \code{tree==TRUE}) `id` (a
##'   unique ID within each simulation for each individual element of the
##'   chain), `ancestor` (the ID of the ancestor of each element) and
##'   `generation`. 
##' @author Sebastian Funk
##' @export
##' @examples
##' chain_sim(n=5, "pois", "size", lambda=0.5)
chain_sim <- function(n, offspring, stat = c("size", "length"), infinite = Inf,
                      tree = FALSE, serial, t0 = 0, tf = Inf, ...) {

    stat <- match.arg(stat)

    ## first, get random function as given by `offspring`
    if (!is.character(offspring)) {
        stop("object passed as 'offspring' is not a character string.")
    }

    roffspring_name <- paste0("r", offspring)
    if (!(exists(roffspring_name)) || !is.function(get(roffspring_name))) {
        stop("Function ", roffspring_name, " does not exist.")
    }

    if (!missing(serial)) {
        if (!is.function(serial)) {
            stop("The `serial` argument must be a function.")
        }
        if (!missing(tree) && tree == FALSE) {
            stop("The `serial` argument can't be used with `tree==FALSE`.")
        }
        tree <- TRUE
    } else if (!missing(tf)) {
        stop("The `tf` argument needs a `serial` argument.")
    }

    stat_track <- rep(1, n) ## track length or size (depending on `stat`)
    n_offspring <- rep(1, n) ## current number of offspring
    sim <- seq_len(n) ## track chains that are still being simulated

    ## initialise data frame to hold the trees
    if (tree) {
        generation <- 1L
        tdf <-
            data.frame(n = seq_len(n),
                       id = 1L,
                       ancestor = NA_integer_,
                       generation = generation)

        ancestor_ids <- rep(1, n)
        if (!missing(serial)) {
            tdf$time <- t0
            times <- tdf$time
        }
    }

    ## next, simulate n chains
    while (length(sim) > 0) {
        ## simulate next generation
        next_gen <- get(roffspring_name)(n=sum(n_offspring[sim]), ...)
        if (any(next_gen %% 1 > 0)) {
            stop("Offspring distribution must return integers")
        }

        ## record indices corresponding the number of offspring
        indices <- rep(sim, n_offspring[sim])

        ## initialise number of offspring
        n_offspring <- rep(0, n)
        ## assign offspring sum to indices still being simulated
        n_offspring[sim] <- tapply(next_gen, indices, sum)

        ## track size/length
        if (stat=="size") {
            stat_track <- stat_track + n_offspring
        } else if (stat=="length") {
            stat_track <- stat_track + pmin(1, n_offspring)
        }

        ## record times/ancestors (if tree==TRUE)
        if (tree && sum(n_offspring[sim]) > 0) {
            ancestors <- rep(ancestor_ids, next_gen)
            current_max_id <- unname(tapply(ancestor_ids, indices, max))
            indices <- rep(sim, n_offspring[sim])
            ids <- rep(current_max_id, n_offspring[sim]) +
                unlist(lapply(n_offspring[sim], seq_len))
            generation <- generation + 1L
            new_df <-
                data.frame(n = indices,
                           id = ids,
                           ancestor = ancestors,
                           generation = generation)
            if (!missing(serial)) {
                times <- rep(times, next_gen) + serial(sum(n_offspring))
                current_min_time <- unname(tapply(times, indices, min))
                new_df$time <- times
            }
            tdf <- rbind(tdf, new_df)
        }

        ## only continue to simulate chains that offspring and aren't of
        ## infinite size/length
        sim <- which(n_offspring > 0 & stat_track < infinite)
        if (length(sim) > 0) {
            if (!missing(serial)) {
                ## only continue to simulate chains that don't go beyond tf
                sim <- intersect(sim, unique(indices)[current_min_time < tf])
            }
            if (tree) {
                if (!missing(serial)) {
                    times <- times[indices %in% sim]
                }
                ancestor_ids <- ids[indices %in% sim]
            }
        }
    }

    if (tree) {
        if (!missing(tf)) {
            tdf <- tdf[tdf$time < tf, ]
        }
        rownames(tdf) <- NULL
        return(tdf)
    } else {
        stat_track[stat_track >= infinite] <- Inf
        return(stat_track)
    }
}

