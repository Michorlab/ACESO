
models_included <- function() {
  message("ACESO internal library of PK models:")
  models <- gsub(".cpp","",dir(file.path(path.package("ACESO"), "PK_models")))
  print(models)
  return(invisible(NULL))
}
##' Internal model library from ACESO
##' 
##' @param list logical:list available models
##' @export
##' 
##' @details
##' Call \code{model_library(list=TRUE)} to list available models.  Once the model 
##' is loaded (see examples below), call \code{(mod)@code} to see
##' model code and equations.
##' 
##' 
##' @examples
##' \dontrun{
##' model_library(list=TRUE)
##' mod <- mread("1cmt_iv", model_library())
##' mod
##' mod <- mread("1cmt_ev", model_library()) 
##' mod@code
##' 
##' } 

model_library <- function(list=FALSE)  {
  if(list) {
    return(models_included())
  }
  return(file.path(path.package("ACESO"), "PK_models"))
}