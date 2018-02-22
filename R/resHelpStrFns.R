# Set up the options for the package:

optEnv <- new.env()
optDefaults <- list(
  'digits'=2,
  'dfDigits'=0,
  'pDigits'=2,
  'confLvl'=.95,
  'asEqn'=T,
  'template'='f_eta2_p'
)
for (o in names(optDefaults)) {
  assign(o, optDefaults[[o]], envir=optEnv)
}

rm(optDefaults, o)

getOpts <- function (arguments=NULL, ...) {
  #' A way to quickly get options for functions.
  #' 
  #' @param arguments A character vector of arguments to access. If NULL,
  #' the function processes all available arguments.
  #' @param ... Arguments that override the defaults. If an argument is provided
  #' via ellipses from a caller function, these are assigned the provided value 
  #' instead of the default value.
  #' 
  #' @details 
  #' If the 
  
  # If there are ellipsis arguments, then assign them first.
  
  overrideArgs <- list(...)
  for (arg in names(overrideArgs)) {
    assign(arg, overrideArgs[[arg]], envir = parent.frame())
  }
  
  if (is.null(arguments)) {
    arguments <- ls(optEnv)
  }
  
  for (arg in arguments) {
    if (!exists(arg, envir = parent.frame())) {
      #browser()
      # If it does not exist, then assign it from the defaults.
      assign(arg, get(arg, envir = optEnv), envir=parent.frame())
    }
  }
}

# Actual Helper Functions

resHelpPrep <- function(param, modObj) {
  #' Converts a model object to a summary object if necessary, validates
  #' the parameter information, and grabs options if necessary.
  #' 
  #' @param param A string containing the parameter of interest. If equal to 
  #' "total" then the overall F is used. The shortcut "int" works as well.
  #' @param modObj Either an \code{lm} object, or a \code{summary.lm} object.
  #' @param digits General number of digits to round to.
  #' @param asEqn Should the statement be formatted as an equation?
  #' @param dfDigits Digits to round degrees of freedom to.
  #' @param pDigits Digits to round p-value to.
  #' 
  
  # Make sure we have at least an lm object or a summary
  if (!class(modObj) %in% c('lm', 'summary.lm')) {
    stop('need an lm object or a summary object')
  }
  
  # make a summary object if necessary
  if (!class(modObj)=='summary.lm') {
    lmSumObj <- summary(modObj)
  } else {
    lmSumObj <- modObj
  }
  rm(modObj)
  
  # validate the param argument
  if (param == 'int') {
    param <- '(Intercept)'
  } else if (!all(param %in% rownames(lmSumObj$coefficients))) {
    stop("param must be in the model")
  }
  
  return(list(param=param, 
              lmSumObj=lmSumObj
              ))
}

fStr <- function (param, modObj, ...) {
  #' Constructs an APA formatted string for an F value.
  #' 
  #' @param param The parameter of interest. Can be "total" for the omnibus
  #' f-value, or "int" as a shortcut for "(Intercept)".
  #' @param modObj Either an \code{lm} or \code{summary.lm} (faster) object.
  #' @param ... Options.
  #' @export fStr
  
  # Get default values we need if not provided
  getOpts(...)
  
  # Get model info
  validated <- resHelpPrep(param, modObj)
  param = validated$param; lmSumObj <- validated$lmSumObj
  
  # Get the right f value and degrees of freedom
  if (param == 'total') {
    if (lmSumObj$df[1] > 1) {
      f <- lmSumObj$fstatistic['value']
      df1 <- lmSumObj$fstatistic['numdf']
    } else {
      f <- (lmSumObj$coefficients['(Intercept)','t value']) ^ 2
      df1 <- 1
    }
  } else {
    f <- (lmSumObj$coefficients[param,'t value']) ^ 2
    df1 <- 1
  }
  df2 <- lmSumObj$df[2]
  
  # round all the numbers
  f <- formatC(f, digits=digits, format='f')
  df1 <- formatC(df1, digits=dfDigits, format='f') 
  df2 <- formatC(df2, digits=dfDigits, format='f')
  
  # make a string!
  str <- paste0('F(', df1, ',', df2, ')=', f)
  if (asEqn) {
    str <- paste0('$', str, '$')
  }
  return(str)
}

eta2Str <- function (param, modObj, ...) {
  #' Constructs an APA formatted string for an \eqn{\eta^2_p} value. Note that
  #' this needs to be an equation for the character to display correctly in 
  #' Markdown (\eqn{LaTeX}).
  #' 
  #' @param param The parameter of interest. Can be "total" for the omnibus
  #' f-value, or "int" as a shortcut for "(Intercept)".
  #' @param modObj Either an \code{lm} or \code{summary.lm} (faster) object.
  #' @param ... Options.
  #' @export fStr
  
  # Get default values we need if not provided
  getOpts(...)
  
  # Get model info
  validated <- resHelpPrep(param, modObj)
  param = validated$param; lmSumObj <- validated$lmSumObj
  
  # Get the right R^2 value
  fval <- (lmSumObj$coefficients[param,'t value']) ^ 2
  rdf <- lmSumObj$df[2]
  r2 <- fval/(fval+rdf)
  
  # round all the numbers
  r2 <- formatC(r2, digits=digits, format='f')
  
  # make a string!
  str <- paste0('\\eta^2_p=', r2)
  if (asEqn) {
    str <- paste0('$', str, '$')
  }
  return(str)
}

pStr <- function (param, modObj, ...) {
  #' Constructs an APA formatted string for a p-value.
  #' 
  #' @param param The parameter of interest. Can be "total" for the omnibus
  #' f-value, or "int" as a shortcut for "(Intercept)".
  #' @param modObj Either an \code{lm} or \code{summary.lm} (faster) object.
  #' @param ... Options.
  #' @export pStr
  
  # Get options.
  getOpts(...)
  
  prep <- resHelpPrep (param, modObj)
  param <- prep$param; lmSumObj <- prep$lmSumObj
  
  # Get the right p-value.
  if (param=='total') {
    pVal <- 
      1 - pf(lmSumObj$fstatistic['value'], 
             lmSumObj$fstatistic['numdf'],
             lmSumObj$fstatistic['dendf'])
  } else {
    pVal <- lmSumObj$coefficients[param, 'Pr(>|t|)']
  }
  
  # Make a string
  if (pVal < .001) {
    baseStr <- "p<.001"
  } else if (pVal < .01 & pDigits <= 2) {
    baseStr <- "p<.01"
  } else {
    baseStr <-  paste0('p=', formatC(pVal, digits=pDigits, format='f'))
  }
  
  if (asEqn) {
    str <- paste0('$', baseStr, '$')
  } else {
    str <- baseStr
  }
  
  return(str)
}

tStr <- function (param, modObj, ...) {
  #' Constructs an APA formatted string for a t-value.
  #' 
  #' @param param The parameter of interest. Can be "int" as a shortcut 
  #' for "(Intercept)".
  #' @param modObj Either an \code{lm} or \code{summary.lm} (faster) object.
  #' @param ... Options.
  #' @export tStr
  
  # Get options
  getOpts(...)
  
  prep <- resHelpPrep (param, modObj)
  param <- prep$param; lmSumObj <- prep$lmSumObj
  
  # Get the right p-value.
  if (param=='total') {
    stop('t-values do not make sense for omnibus tests! Use F instead.')
  } else {
    tVal <- lmSumObj$coefficients[param, 't value']
    df <- lmSumObj$df[2]
  }
  
  # Make a string
  baseStr <- paste0('t(', 
                    formatC(df, digits=dfDigits, format='f'), 
                    ')=', 
                    formatC(tVal, digits=digits, format='f')
  )
  
  if (asEqn) {
    str <- paste0('$', baseStr, '$')
  } else {
    str <- baseStr
  }
  
  return(str)
}

mStr <- function (v, na.rm=F, digits=NULL, asEqn=NULL) {
  #' Prints an APA formatted mean.
  #' 
  #' @param v Either a mean value or a vector. If a vector, then \code{mean()} 
  #' is called on it.
  #' @param na.rm If v is a vector, this is passed to \code{mean()}. 
  #' Default: FALSE.
  #' @export mStr 
  
  # Get options
  getOpts()
  
  if (length(v) > 1) {
    mn <- mean(v, na.rm=na.rm)
  } else {
    mn <- v
  }
  
  # Make a string!
  str <- paste0('M=', formatC(mn, digits=digits), format='f')
  
  if (asEqn) {
    str <- paste0('$', str, '$')
  }
  
  return(str)
}

sdStr <- function (v, na.rm=F, digits=NULL, asEqn=NULL) {
  #' Prints an APA formatted sd.
  #' 
  #' @param v Either a mean value or a vector. If a vector, then \code{sd()} 
  #' is called on it.
  #' @param na.rm If v is a vector, this is passed to \code{sd()}. 
  #' Default: FALSE.
  #' @export sdStr 
  
  if (length(v) > 1) {
    sd <- sd(v, na.rm=na.rm)
  } else {
    sd <- v
  }
  
  # Get options
  getOpts()
  
  # Make a string!
  str <- paste0('SD=', formatC(sd, digits=digits, format='f'))
  
  if (asEqn) {
    str <- paste0('$', str, '$')
  }
  
  return(str)
}

confIntSum <- function (object, parm, confLvl=NULL) {
  #' Confidence interval for summary objects. Derived directly from 
  #' \code{confint.lm()}.
  #' 
  #' @param object Summary object of lm model.
  #' @param parm a specification of which parameters are to be given confidence 
  #' intervals, either a vector of numbers or a vector of names. If missing, 
  #' all parameters are considered.
  #' @param confLvl Alpha value.
  #' 
  
  getOpts()
  
  cf <- coef(object)[,"Estimate"]
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - confLvl)/2
  a <- c(a, 1 - a)
  fac <- qt(a, object$df[2]) # residual degrees of freedom for the model.
  #pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             a))
  ses <- sqrt(diag(vcov(object)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}


confIntStr <- function (param, modObj, digits=NULL, confLvl=NULL, asEqn=NULL) {
  #' Constructs an APA formatted string for a confidence interval around
  #' a parameter estimate from a model.
  #' 
  #' @param param The parameter of interest. Can be "int" as a shortcut 
  #' for "(Intercept)".
  #' @param modObj Either an \code{lm} or \code{summary.lm} (faster) object.
  #' @param digits Number of digits to round to.
  #' @param confLvl The alpha for confidence interval estimation.
  #' @param asEqn Should the string be an equation?
  #' @export confIntStr
  #' 
  
  getOpts()
  
  prep <- resHelpPrep(param, modObj)
  param <- prep$param; sumObj <- prep$lmSumObj
  
  interval <- confIntSum(sumObj, param)
  str <- paste0('95% CI [', 
                formatC(interval[1], digits=digits, format='f'),
                ', ',
                formatC(interval[2], digits=digits, format='f'),
                ']'
                )
  if (asEqn) {
    str <- paste0('$', str, '$')
  }
  return(str)
}


fpStr <- function (param, modObj, asEqn=NULL, digits=NULL, dfDigits=NULL) {
  #' Prints an APA formated string of the form:
  #' \code{F(df1,df2),p[=,<]pval}.
  #' 
  #' @param param The parameter of interest. The shortcut "int" is recognized
  #' as well as the "total" term.
  #' @param modObj The \code{lm} or \code{summary.lm} object.
  #' @param asEqn Should output be a markdown equation?
  #' @param ... Other options.
  #' 
  #' @export fpStr
  
  # Get options
  getOpts()
  
  fstr <- fStr(param, modObj, asEqn=F, digits=digits, dfDigits=dfDigits)
  pstr <- pStr(param, modObj, asEqn=F, pDigits)
  
  if (asEqn) {
    str <- paste0('$', fstr, '$, $', pstr, '$')
  } else {
    str <- paste0(fstr, ', ', pstr)
  }
  
  return(str)
}

tpStr <- function (param, modObj, asEqn=NULL, ...) {
  #' Prints an APA formated string of the form:
  #' \code{t(df),p[=,<]pval}.
  #' 
  #' @param param The parameter of interest. The shortcut "int" is recognized.
  #' @param modObj The \code{lm} or \code{summary.lm} object.
  #' @param asEqn Should output be a markdown equation?
  #' @param ... Other options.
  #' 
  #' @export tpStr
  
  getOpts()
  
  tstr <- tStr(param, modObj, asEqn=F, digits=digits, dfdigits=dfdigits)
  pstr <- pStr(param, modObj, asEqn=F, pdigits=pdigits)
  
  if (asEqn) {
    str <- paste0('$', tstr, '$, $', pstr, '$')
  } else {
    str <- paste0(tstr, ', ', pstr)
  }
  
  return(str)
}

printStats <- function (param, modObj, ...) {
  #' @title Print a Formatted Stats String
  #' 
  #' @param param The parameter of interest. The shortcut "int" is recognized.
  #' @param modObj The \code{lm} or \code{summary.lm} object.
  #' @param template A string specifying which stats to include in this string. 
  #' See details for how to construct a string. If left as \code{NULL}, the 
  #' default template is used.
  #' @param ... Other options.
  #' 
  #' @details 
  #' The \code{template} argument is a string containing a sequence of the
  #' following characters, each separated by an underscore (always lowercase):
  #' \itemize{
  #'   \item f: An 1-degree-of-freedom F-statistic, e.g., \eqn{F(1,36)=4.2}
  #'   \item t: A t-statistic, e.g., \eqn{t(38)=2.1}
  #'   \item p: A p-value, e.g., \eqn{p=.03} or \eqn{p<.01} depending on the
  #'   value
  #'   \item eta2: An \eqn{\eta^2_p} value, e.g., \eqn{\eta^2_p=0.056}
  #'   \item m: A mean, e.g., \eqn{M=.45}.
  #'   \item sd: A standard deviation, e.g., \eqn{SD=.34}
  #'   \item ci: A confidence interval, e.g., \eqn{95%CI[.01, 2.34]}
  #' }
  #' Note that the same set of \code{...} arguments is passed to every function 
  #' call. Otherwise, the strings are built with default argument values as
  #' normal.
  #' 
  #' @export printStats
  
  getOpts(...)
  
  fnNms <- list(
    f=fStr,
    t=tStr,
    p=pStr,
    ci=confIntStr,
    m=mStr,
    sd=sdStr,
    eta2=eta2Str
  )
  fnStr <- strsplit(template, '_')[[1]]
  strLst <- sapply(fnStr, function (fnNm) {
    fnNms[[fnNm]](param = param, modObj = modObj, asEqn=F, ...)
  })
  
  if (asEqn) {
    strLst <- sapply(strLst, function (str) {
      return(paste0('$', str, '$'))
    })
  }
  
  mainStr <- Reduce(function (str1, str2) {
    return(paste0(str1, ', ', str2))
  }, strLst)
  
  return(mainStr)
}

