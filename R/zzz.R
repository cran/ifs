if(R.Version()$major > 1) { require(stats) 
}else{  require(stepfun) }

.First.lib <- function(lib, pkg) library.dynam("ifs", pkg, lib) 
