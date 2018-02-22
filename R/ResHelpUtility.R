#' Save All Model Objects to a File
#' 
#' This function is intended to save all objects of class "lm" to a specified
#' file. Note that by changing the \code{objClass} argument, it may be used to
#' save any kind of object to a file. You may also specify the environment in 
#' which to search for objects.
#'
#' @param file A string representing a filename to pass to \code{save()}.
#' @param envir Any way of specifying an environment in \code{ls()} and
#' \code{get()}.
#' @param objClass A character string of the class of object to be selected.
#' @param ... Other options to pass to \code{save}
#'
#' @export saveAllModelObjects
#'
saveAllModelObjects <- function (file, envir=".GlobalEnv", objClass="lm", ...) {
  if (missing(file)) {
    stop("No filename provided")
  }
  envir <- as.environment(envir)
  allObjNms <- ls(name=envir)
  modObjNms <- allObjNms[sapply(allObjNms, 
                                function(obj) {
                                  curClass <- class(get(obj, envir = envir))
                                  return(objClass %in% c(curClass))
                                  })]
  
  save(list = modObjNms, 
       envir = envir, 
       file = file)
}