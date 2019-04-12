### Define selfStart functions for fitting reaction rate data
### as a function of substrate and inhibitor concentrations to
### models for competitive, uncompetitive and mixed inhibition.

#' Model different types of enzyme kinetic behaviors for nonlinear
#'   optimization of rate data.
#'
#' The selfStart functions provided by the package are
#' \itemize{
#'   \item \code{\link{SScompinh}} Model for competitive inhibition.
#'   \item \code{\link{SSuncompinh}} Model for uncompetitive inhibition.
#'   \item \code{\link{SSmixedinh}} Model for mixed inhibition.
#' }
#'
#' The functions are intended to be used as the right-hand side of
#' formulas in \code{\link[stats]{nls}} or similar functions. See the
#' documentation of the functions for examples.
"_PACKAGE"

#' Supplies a selfStart function for modeling competitive inhibition of
#'   enzyme kinetics.
#'
#' @param sub A vector of substrate concentrations.
#' @param inh A vector of inhibitor concentrations.
#' @param Vmax The maximum initial reaction rate.
#' @param Km The Michaelis constant for the enzyme-substrate complex.
#' @param Kc The dissociation constant for the enzyme-inhibitor complex.
#' @return A \code{\link{selfStart}} function for use in nonlinear
#'           optimization.
#' @seealso \code{\link{selfStart}}, \code{\link[stats]{nls}}
#' @export
#' @examples
#' \dontrun{compnls <- nls(Rate~SScompinh(Sub.Con,Inh.Con,Vm,Km,Kcomp),
#'                         data=Enzymeframe)}
SScompinh <- selfStart(~ Vmax*sub/(Km*(1+inh/Kc)+sub),
                       function(mCall, data, LHS)
                       {
                           ## Calculate vectors for linear model
                           invrate <- 1/eval(LHS, envir=data)
                           invsub <- 1/eval(mCall[["sub"]], envir=data)
                           inhsub <- eval(mCall[["inh"]], envir=data)/eval(mCall[["sub"]], envir=data)
                           ## Calculate the coefficients of the linear
                           ## model directly from the linear model.
                           inhcoef <- coef(lm(invrate~invsub+inhsub))
                           pars <- c(1/inhcoef[["(Intercept)"]],
                                     inhcoef[["invsub"]]/inhcoef[["(Intercept)"]],
                                     inhcoef[["invsub"]]/inhcoef[["inhsub"]])
                           setNames(pars, mCall[c("Vmax","Km","Kc")])
                       }, c("Vmax","Km","Kc"))

#' Supplies a selfStart function for modeling uncompetitive inhibition of
#'   enzyme kinetics.
#'
#' @param sub A vector of substrate concentrations.
#' @param inh A vector of inhibitor concentrations.
#' @param Vmax The maximum initial reaction rate.
#' @param Km The Michaelis constant for the enzyme-substrate complex.
#' @param Ku The dissociation constant for the ternary complex.
#' @return A \code{\link{selfStart}} function for use in nonlinear
#'           optimization.
#' @seealso \code{\link[stats]{selfStart}}, \code{\link[stats]{nls}}.
#' @export
#' @examples
#' \dontrun{x <- nls(Rate~SSuncompinh(Sub.Con,Inh.Con,Vm,Km,Ku),
#'                   data=EnzymeFrame)}
SSuncompinh <- selfStart(~ Vmax*sub/(Km+(1+inh/Ku)*sub),
                       function(mCall, data, LHS)
                       {
                           ## Calculate vectors for linear model
                           invrate <- 1/eval(LHS, envir=data)
                           invsub <- 1/eval(mCall[["sub"]], envir=data)
                           inhibitor <- eval(mCall[["inh"]], envir=data)
                           ## Calculate coefficients of linear model
                           ## directly.
                           inhcoef <- coef(lm(invrate~invsub+inhibitor))
                           pars <- c(1/inhcoef[["(Intercept)"]],
                                     inhcoef[["invsub"]]/inhcoef[["(Intercept)"]],
                                     inhcoef[["(Intercept)"]]/inhcoef[["inhibitor"]])
                           setNames(pars, mCall[c("Vmax","Km","Ku")])
                       }, c("Vmax","Km","Ku"))

#' Supplies a selfStart function to model mixed inhibition of enzyme
#'   kinetics.
#'
#' @param sub A vector of substrate concentrations.
#' @param inh A vector of inhibitor concentrations.
#' @param Vmax The maximum initial rate of reaction.
#' @param Km The Michaelis constant for the enzyme-substrate complex.
#' @param Kc The dissociation constant for the enzyme-inhibitor complex.
#' @param Ku The dissociation constant for the ternary complex.
#' @return A \code{\link{selfStart}} function for use in nonlinear
#'          optimization.
#' @seealso \code{\link[stats]{selfStart}}, \code{\link[stats]{nls}}
#' @export
#' @examples
#' \dontrun{mixed <- nls(Rate~SSmixedinh(Sub.Con,Inh.Con,Vm,Km,Kc,Ku),
#'                       data=EnzymeFrame)}
SSmixedinh <- selfStart(~ Vmax*sub/(Km*(1+inh/Kc)+sub*(1+inh/Ku)),
                       function(mCall, data, LHS)
                       {
                           ## Calculate vectors for linear model
                           invrate <- 1/eval(LHS, envir=data)
                           invsub <- 1/eval(mCall[["sub"]], envir=data)
                           inhsub <- eval(mCall[["inh"]], envir=data)/eval(mCall[["sub"]], envir=data)
                           inhibitor <- eval(mCall[["inh"]], envir=data)
                           ## Calculate coefficients of linear model
                           ## directly.
                           inhcoef <- coef(lm(invrate~invsub+inhsub+inhibitor))
                           pars <- c(1/inhcoef[["(Intercept)"]],
                                     inhcoef[["invsub"]]/inhcoef[["(Intercept)"]],
                                     inhcoef[["invsub"]]/inhcoef[["inhsub"]],
                                     inhcoef[["(Intercept)"]]/inhcoef[["inhibitor"]])
                           setNames(pars, mCall[c("Vmax","Km","Kc","Ku")])
                       }, c("Vmax","Km","Kc","Ku"))
