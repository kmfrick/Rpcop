pcop <- function(x,Ch=1.5,Cd=.3,plot.true=FALSE,...){
	proyec <- as.data.frame(pcop_backend(x, Cd, Ch));

	pr.curve <- proyec[,-1]

		p <- dim(x)[2]

		names(pr.curve) <- c("param", "dens", "span", "orth.var", 
				paste("pop", 1:p, sep = ""), 
				paste("pr.dir", 1:p, sep = ""))
		# pr.curve: principal curve of oriented points 
		#           with original format (format 1)

		# pc: principal curve of oriented points 
		#     with format as in HS's "princurv" library (format 2)
		s<-as.matrix(pr.curve[,5:(5+p-1)])
		#pc <- get.lam(x,s,...) # get.lam is deprecated, please use project_to_curve instead
		pc <- project_to_curve(x,s)

		if (plot.true) {
			adjust.range <- function(x,cte=1){
				rgx<-range(x)
					extra <- diff(rgx)*(cte-1)/2
					argx1<-rgx[1]-extra
					argx2<-rgx[2]+extra
					return(c(argx1,argx2))
			}
			plot(x[, 1:2], xlim = adjust.range(x[, 1], 1.1), ylim = adjust.range(x[,2], 1.1), col=8)
				lines(pc$s[pc$ord, 1:2],...)
		}


	return(list(pcop.f1=pr.curve,pcop.f2=pc))
}

