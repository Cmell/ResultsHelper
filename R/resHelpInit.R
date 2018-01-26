# Set up the options for the package:

optEnv <- new.env()
optDefaults <- list(
  'resHelp-digits'=2,
  'resHelp-dfDigits'=0,
  'resHelp-pDigits'=2,
  'resHelp-asEqn'=T
  )
for (o in names(optDefaults)) {
  assign(o, optDefaults[[o]], envir=optEnv)
}
rm(optDefaults, o)

getOpts <- function (valLst, ...) {
  retVal <- lapplay(valLst, function (v) {
    if (missing(v)) {get(v, envir=optEnv)}
  })
  return(retVal)
}