# This file contains R functions that are used in the simulation study.

# This function prints the proportions of CI contenting the 
# true value for a vector of parameters
f.coverage.vec <- function(out, theta0){
  out      <- out[!sapply(out, is.null)]
  rbind("90%.normal1" = apply(sapply(out, function(x) theta0 >= x$quant1[,"5%"] & theta0 <= x$quant1[,"95%"]), 1, mean),
        "95%.normal1" = apply(sapply(out, function(x) theta0 >= x$quant1[,"2.5%"] & theta0 <= x$quant1[,"97.5%"]), 1, mean),
        "99%.normal1" = apply(sapply(out, function(x) theta0 >= x$quant1[,"0.5%"] & theta0 <= x$quant1[,"99.5%"]), 1, mean),
        "90%.normal2" = apply(sapply(out, function(x) theta0 >= x$quant2[,"5%"] & theta0 <= x$quant2[,"95%"]), 1, mean),
        "95%.normal2" = apply(sapply(out, function(x) theta0 >= x$quant2[,"2.5%"] & theta0 <= x$quant2[,"97.5%"]), 1, mean),
        "99%.normal2" = apply(sapply(out, function(x) theta0 >= x$quant2[,"0.5%"] & theta0 <= x$quant2[,"99.5%"]), 1, mean),
        "90%.approxi" = apply(sapply(out, function(x) theta0 >= x$quant[,"95%"] & theta0 <= x$quant[,"5%"]), 1, mean),
        "95%.approxi" = apply(sapply(out, function(x) theta0 >= x$quant[,"97.5%"] & theta0 <= x$quant[,"2.5%"]), 1, mean),
        "99%.approxi" = apply(sapply(out, function(x) theta0 >= x$quant[,"99.5%"] & theta0 <= x$quant[,"0.5%"]), 1, mean))
}

# This function prints the proportions of CI contenting the 
# true value for a scalar
f.coverage.sca <- function(out, theta0){
  out      <- out[!sapply(out, is.null)]
  rbind("90%.normal1" = mean(sapply(out, function(x) theta0 >= x$quant1["5%"] & theta0 <= x$quant1["95%"])),
        "95%.normal1" = mean(sapply(out, function(x) theta0 >= x$quant1["2.5%"] & theta0 <= x$quant1["97.5%"])),
        "99%.normal1" = mean(sapply(out, function(x) theta0 >= x$quant1["0.5%"] & theta0 <= x$quant1["99.5%"])),
        "90%.normal2" = mean(sapply(out, function(x) theta0 >= x$quant2["5%"] & theta0 <= x$quant2["95%"])),
        "95%.normal2" = mean(sapply(out, function(x) theta0 >= x$quant2["2.5%"] & theta0 <= x$quant2["97.5%"])),
        "99%.normal2" = mean(sapply(out, function(x) theta0 >= x$quant2["0.5%"] & theta0 <= x$quant2["99.5%"])),
        "90%.approxi" = mean(sapply(out, function(x) theta0 >= x$quant["95%"] & theta0 <= x$quant["5%"])),
        "95%.approxi" = mean(sapply(out, function(x) theta0 >= x$quant["97.5%"] & theta0 <= x$quant["2.5%"])),
        "99%.approxi" = mean(sapply(out, function(x) theta0 >= x$quant["99.5%"] & theta0 <= x$quant["0.5%"])))
}


