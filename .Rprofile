args = commandArgs(trailingOnly = T)

head <- function(x, y = 5) { base::print(utils::head(x, y)) }
say <- function(...) { base::print(paste(...)) }

tryCatch({options(width = as.integer(Sys.getenv("COLUMNS")))}, error = function(err) {options(width=236)})

.Last <- function(){
    if (!any(commandArgs() == '--no-readline') && interactive()) {
    	require(utils)
    	try(savehistory(".Rhistory"))
    }
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
