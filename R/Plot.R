Plot.rgasp <- function(object,...) {
  leave_one_out=leave_one_out_rgasp(object)
  
  mean=leave_one_out$mean
  sd=  leave_one_out$sd
  
  resid <- (object@output-mean)/sd
  
  xmin <- min(min(mean), min(object@output))
  xmax <- max(max(mean), max(object@output))
  
  par(mfrow=c(3,1))
  graphics::plot(x = as.numeric(object@output), y = mean,
       xlim = c(xmin, xmax), ylim = c(xmin, xmax),
       xlab = "Exact outputs", ylab = "Fitted outputs",
       main = "Leave-one-out")
  lines(x = c(xmin, xmax), y = c(xmin, xmax))
  graphics::plot(resid, xlab = "Num", ylab = "Standardized residuals",
       main = "Standardized residuals")
  qqnorm(resid, main = "Normal QQ-plot of standardized residuals") 
  qqline(resid)
  par(mfrow = c(1, 1))
  invisible(leave_one_out)
  
}
