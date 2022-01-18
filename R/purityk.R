# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.
# NOTE: The function in this file comes from the semvar Package
# and it was written by Dirk Speelman. In time, I might adapt it to "my style"
# (so I understand it better).
# Code + docs added to semcloud with his permission on 18-jan-2022
# Original version also has a bunch of methods

#' Proportion of same-class nearest neighbours
#'
#' This function takes a square matrix \code{dmx} that contains item by item
#'   distances, and a factor \code{classes}
#'   (with as many items as there are rows, and thus columns,
#'   in \code{dmx}) that assigns a class to each item.
#'   The function returns a measure \bold{q} of how well the distances
#'   in \code{dmx} 'capture' the classification in \code{classes},
#'   where distances are taken to  'capture a classification' to the extent
#'   that items are (immediately) surrounded by other items from the same class,
#'   and not by items from some other class.
#'   Next to an \emph{overall cluster quality} for all the data taken together,
#'   the function also returns the \emph{cluster quality of individual points}
#'   and the  \emph{cluster quality of individual classes}
#'   (as well as the \emph{mean cluster quality over classes}).
#'   All these measures are called \bold{q} in the output of the function.
#'
#'   The \bold{q} measures are calculated as follows:
#'   first, for each item an item-specific cluster quality is calculated.
#'   This is done by calculating the proportion of 'same class items'
#'   among its k nearest neighbours. The higher the measure, the better
#'   the cluster quality for that item. However, what is calculated is not
#'   simply the proportion, but rather the weighted mean of the values
#'   of the k nearest neighbours, where a 'same class neighbour' has value one,
#'   a 'different class neighbour' has value zero, and the weights of the
#'   neighbours can have different settings (see below).
#'   In the default settings, weights decrease linearly with their rank of
#'   'distance from the item', and all weights add up to one. For instance,
#'   if k is one then the weight is 1. If k is 2, then the weights, starting
#'   from the closest nearest neighbour, are .67 and .33. If k is 3, then the
#'   weights are .5, .33, and .17. If k is 4, they are .4, .3, .2, and .1. Etc.
#'
#'   The overall cluster quality of the data is then calculated as the
#'   mean cluster quality of all items. Additionally, the cluster quality
#'   for every class in \code{classes} is calculated as the mean cluster
#'   quality of the items belonging to that class. The mean class quality,
#'   finally, is the mean of all class-specific class quality measures.
#'
#' @param dmx A square matrix containing item by item distances
#' @param classes A factor of the same length as the number of rows and columns
#'   in \code{dmx}; the class in position \emph{i} in \code{classes} is the
#'   class assigned to the item of row \emph{i} and column \emph{i} in \code{dmx}
#' @param k The value of \code{k} that is to be used to identify the \code{k}
#'   nearest neighbours. If \code{k} is not specified, then \code{k} is taken
#'   to be either the total number of items divided by ten (if the number of
#'   items divided by ten is smaller than the size of the smallest class),
#'   or the size of the smallest class minus one (if the size of the smallest
#'   class minus one is smaller than the total number of items divided by ten).
#'   This default behaviour obviously is but a very crude attempt at guessing
#'   a sensible value for \code{k}. Most of the time you probably want to
#'   overrule this default behaviour. If you explicitly specify \code{k},
#'   all value from one up to the total number of items minus one are allowed.
#' @param weights The \code{weights} argument determines how exactly the cluster
#'   quality of a point is derived from the class membership of its k nearest neighbours.
#'   This cluster quality is 'the weighted mean of class membership
#'   values of these neighbours (1=same class as target item; 0=different class)',
#'   with the weights being determined by the \code{weights} argument.
#'   The weights are k numbers, the first of which indicates the weight if
#'   the closest neighbour, the second of which indicates the weight of the
#'   second closest neighbour, etc. The sum of the weights always is one.
#'   When \code{weights} is \code{"linear"}, which is the default situation,
#'   weights decrease linearly as one progresses through the set of neighbours
#'   (starting from the one that is closest to the target item).
#'   When \code{weights} is \code{"s-curve"}, weights decrease as one
#'   progresses through the set of neigbours (starting from the one that is
#'   closest to the target item) according to the s-shape of
#'   \code{y<-(40:-40)/10; plot(1:81, exp(y) / (1 + exp(y)), type="l")},
#'   but with the actual weights rescaled so that they add up to one.
#'   Finally, when \code{weight} is \code{"none"}, all connections in the path
#'   have equal weight. The actual weights that are used in a call to
#'   \code{separationkNN()} can be found in the \code{weights} components in its output.
#'
#' @return An object of the class \code{clustqualkNN},
#'   which is a list containing at least the following components:
#'   \item{globqual}{The \emph{global cluster quality} \bold{q}}
#'   \item{meanclassqual}{The mean of all \emph{class-specific cluster quality values} \bold{q}}
#'   \item{classqual}{A table with for each class its \emph{class-specific clusters quality} \bold{q}}
#'   \item{pointqual}{A numeric vector with for each item its \emph{item-specific cluster quality} \bold{q}}
#'   \item{weights}{A numeric vector with the weights that were used}
#'   \item{k}{A number indication which \code{k} was used}
#'
#' @export
#'
#' @examples
#' # we create a 'point cloud', with points belonging to two classes
#' points <- rbind(matrix(rnorm(100, 2, 2), ncol=2),
#'                 matrix(rnorm(100, 4, 2), ncol=2))
#'                 dst <- dist(points, diag=TRUE, upper=TRUE)
#'                 classes <- as.factor(rep(c("a","b"), c(50, 50)))
#' # we analyse the cluster quality, letting the procedure choose k
#' fitkNN <- separationkNN(dst, classes)
#' summary(fitkNN)
#' fitkNN$globqual        # global cluster quality
#' fitkNN$meanclassqual   # mean class quality
#' fitkNN$classqual       # class-specific quality
#'
#' # we analyse the cluster quality, setting k to 25
#' fitkNN <- separationkNN(dst, classes, k=25)
#' summary(fitkNN)
separationkNN <- function (dmx, classes, k = NULL, weights = c("linear",
                                                         "s-curve", "none"))
{
  if (!is.factor(classes)) {
    classes <- as.factor(classes)
  }

  classlevels <- levels(classes)
  classfreqs <- table(classes)
  sscm1 <- min(classfreqs) - 1

  if (sscm1 < 1) {
    stop("Error: all classes must have at least two members.")
  }

  if (length(classlevels) < 2) {
    stop("Error: there must be at least two classes.")
  }

  if (!is.matrix(dmx)) {
    dmx <- as.matrix(dmx)
  }

  nr <- nrow(dmx)
  diag(dmx) <- NA

  if (is.null(k)) {
    k <- min(round(nr/10), sscm1)
  }
  k <- min(max(k, 1), nr - 1)

  if (is.null(weights)) {
    weights <- "linear"
  } else {
    weights <- weights[[1]]
  }

  if (k == 1) {
    weights <- 1
  } else {
    if (weights == "linear") {
      weights <- k:1/(k * (k + 1)/2)
    } else if (weights == "s-curve") {
      weights <- 8 * ((k - 1):0)/(k - 1) - 4
      weights <- exp(weights)/(1 + exp(weights))
      weights <- weights/sum(weights)
    } else {
      weights <- rep(1/k, k)
    }
  }

  retval <- list()
  nn <- t(apply(dmx, 1, order))
  nn.vect <- as.vector(nn)
  nn.class <- matrix(classes[nn.vect], nrow = nr)
  row.class <- matrix(rep(classes), rep(nr, nr), nrow = nr)
  same.class <- nn.class == row.class
  if (k == 1) {
    pointqual <- as.numeric(same.class[, 1])
  } else {
    pointqual <- apply(same.class[, 1:k], 1, stats::weighted.mean,
                       weights)
  }

  globqual <- mean(pointqual)
  retval[["classfreqs"]] <- classfreqs
  retval[["classqual"]] <- rep(NA, length(classlevels))
  names(retval[["classqual"]]) <- classlevels
  retval[["pointclass"]] <- classes
  retval[["pointqual"]] <- pointqual
  names(retval[["pointqual"]]) <- rownames(dmx)
  for (i in 1:length(classlevels)) {
    level <- classlevels[i]
    subpointqual <- pointqual[classes == level]
    retval[["classqual"]][i] <- mean(subpointqual)
  }
  retval[["globqual"]] <- mean(pointqual)
  retval[["meanclassqual"]] <- mean(retval[["classqual"]])
  retval[["k"]] <- k
  retval[["weights"]] <- weights
  # class(retval) <- "clustqualkNN"
  return(retval)
}
