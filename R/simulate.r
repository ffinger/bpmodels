##' Simulate chains using a branching process
##'
##' @param n number of simulations to run.
##' @param offspring offspring distribution, given as the function used to
##'     generate the number of offspring in each generation, e.g. `rpois` for
##'     Poisson distributed offspring
##' @param stat statistic to calculate ("size" or "length" of chains)
##' @param infinite a size or length from which the size/length is to be
##'     considered infinite
##' @param tree return the tree of infectors
##' @param ... parameters of the offspring distribution
##' @return a vector of sizes/lengths (if \code{tree==FALSE}), or a list of trees (if \code{{tree==TRUE}})
##' @author Sebastian Funk
##' @export
##' @examples
##' chain_sim(n=5, rpois, "size", lambda=0.5)
chain_sim <- function(n, offspring, stat = c("size", "length"), infinite = Inf,
                      tree=FALSE, ...) {

    stat <- match.arg(stat)

    ## first, get random function as given by `offspring`
    if (!is.function(offspring)) {
        stop("object passed as 'offspring' is not a function.")
    }

    ## next, simulate n chains
    dist <- c()
    if (tree) trees <- list()

    ## run simulations
    for (i in seq_len(n)) {
        if (tree) {
            trees[[i]] <- data.frame(id=1L, ancestor=NA_integer_, generation=1L)
        }
        stat_track <- 1 ## track length or size (depending on `stat`)
        state <- 1
        total_size <- 1
        while (state > 0 && stat_track < infinite) {
            new_offspring <- offspring(n=state, ...)
            if (any(new_offspring %% 1 > 0)) {
                stop("Offspring distribution must return integers")
            }
            n_offspring <- sum(new_offspring)
            if (tree && n_offspring > 0) {
                current_gen <- max(trees[[i]]$generation)
                ancestors <- trees[[i]][trees[[i]]$generation == current_gen, "id"]
                new_data <-
                    data.frame(id=max(trees[[i]]$id) + seq_len(n_offspring),
                               ancestor=rep(unlist(ancestors), new_offspring),
                               generation=current_gen + 1)
                trees[[i]] <- rbind(trees[[i]], new_data)
            }
            if (stat=="size") {
                stat_track <- stat_track + n_offspring
            } else if (stat=="length") {
                if (n_offspring > 0) stat_track <- stat_track + 1
            }
            state <- n_offspring
        }
        if (stat_track >= infinite) stat_track <- Inf
        dist[i] <- stat_track
    }

    if (tree) {
        retval <- trees
    } else {
        retval <- dist
    }

    return(retval)
}

