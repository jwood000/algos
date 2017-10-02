#' @title Sieve for Euler's Phi Function
#'
#' @description Quickly generates a list of \code{n} integers such that each index \code{k}, where \code{1 <= k <= n}, represents Euler's Phi function acting on \code{k}.
#' @param n Integer
#' @details Euler's Phi function, or totient function, returns the number of integers less than \code{n} that are relatively prime (coprime) to \code{n} (E.g. \code{phi(10) = 4} as 1, 3, 7, & 9 are coprime to 10).  The formula for an individual number \code{n} is given by
#'  \code{phi(n) = n * (1 - 1/p1) * ... * (1 - 1/pn)}
#' where \code{p1,...,pn} are the prime numbers that divide \code{n}.
#' @return Integer vector of length \code{n}, such that each
#' @aliases
#' @examples
eulerPhiSieve <- function(n = 100L) {
    EulerPhiSieveRcpp(n)
}
