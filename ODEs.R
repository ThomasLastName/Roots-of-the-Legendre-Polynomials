

#
##
###	the increment function in classical fourth order Runge-Kutta
##
#

Increment.RK4 <- function(from, to, f, y.at.from) {

	a <- from
	b <- to
	yhat <- y.at.from
	h <- b-a

	K1 <- h*f(a, yhat)
	K2 <- h*f(a + 0.5*h, yhat + 0.5*K1)
	K3 <- h*f(a + 0.5*h, yhat + 0.5*K2)
	K4 <- h*f(b, yhat + K3)
	
	yhat + (K1 + 2*K2 + 2*K3 + K4)/6

}


#
##
###
####	classical Runge-Kutta of order four
###
##
#

Runge.Kutta <- function(diff.eq, interval, initial.value, step.size=NULL, custom.points=NULL, plot=T, print=NULL, store=NULL, points=FALSE, add=FALSE, calc.dy=FALSE) {
#
#	warnings
#
Custom.Points <- custom.points
dim <- length(initial.value)
if (dim>2 & plot) warning('plot=T Will Be Ignored Because dim(initial.value)>2.')
if (!is.function(diff.eq)) stop('f Must Be A Function Of Two Variables, Equal To The Derivative Of Its Second Argument, i.e., f(t,y)=dy/dt.')
if (!(is.character(store) | is.null(store))) stop('store Argument Must Be a Character String or NULL.')
if (length(interval)<2) stop('Please Make interval of the Form c(a,b) where a<b.')
if (!interval[2]>interval[1]) stop('Please Make interval of the Form c(a,b) where a<b.')
if (!is.numeric(interval)) stop('Please Make interval of the Form c(a,b) where a<b.')
if (length(interval)>2) warning('interval Has Length Greater Than Two; Places After Two Ignored.')
if (!is.logical(plot)) stop('plot Must Be Logical.')
if (!is.logical(add)) stop('add Must Be Logical.')
if (!is.logical(calc.dy)) stop('calc.dy Must Be Logical.')
if (!(is.logical(print) | is.null(print))) stop('print Should Be Either Logical Or NULL.')
test <- diff.eq(interval[1], initial.value)
if (!length(test)==dim) stop('The Dimension Of The Codomain Of diff.eq Must Match The Dimension Of initial.value.')
if (!is.null(dim(test))) stop('We Have !is.null(dim(test)) Where test <- diff.eq(interval[1], initial.value).')
if (!is.numeric(test) & !is.complex(test)) stop('We Have !is.numeric(test) & !is.complex(test) Where test <- diff.eq(interval[1], initial.value).')
if (plot & !calc.dy & dim==1) {
	calc.dy <- TRUE
	warning('calc.dy Set To TRUE Automatically In Order To Construc The Hermite Spline.')
}
#
#	the algorithm
#

	f <- diff.eq
	t <- interval[1]
	y <- initial.value
	dy <- test	
	if (is.null(Custom.Points)) {custom <- FALSE} else {custom <- TRUE}
	if (dim>1) {
		inc <- y
		incinc <- dy
		y <- list(inc[1])
		dy <- list(incinc[1])
		for (i in 1:dim) {
			y[[i]] <- inc[i]
			dy[[i]] <- incinc[i]
			names(y[[i]]) <- paste('y', i, sep='')
			names(dy[[i]]) <- paste('dy', i, sep='')
		}
	}

if (custom) {

	if (!is.null(step.size)) warning('Argument step.size Ignored Because custom.points Was Supplied.')
	t <- Custom.Points

} else {

	if (is.null(step.size)) step.size <- 0.01
	if (!is.null(Custom.Points)) warning('Argument xustom.points Ignored Because No Points Were Supplied.')
	h <- step.size
	t <- seq(interval[1], interval[2], h)
#	if (!(interval[2] %in% t)) t[length(t)+1] <- interval[2]

}
	
	if (dim>1) {
		for (i in 2:length(t)) {
			inc <- Increment.RK4(t[i-1], t[i], f, inc)
			for (j in 1:dim) y[[j]][i] <- inc[j]
			if (calc.dy) {
				dee <- f(t, inc)
				for (j in 1:dim) dy[[j]][i] <- dee[j]
			}
		}
	} else {
		for (i in 2:length(t)) y[i] <- Increment.RK4(t[i-1], t[i], f, y[i-1])
		if (calc.dy) dy <- f(t, y)
	}

#
#	what to do with the results
#
	if (
		plot & ( (dim==1 & any(is.complex(y))) | (dim>1 & any(is.complex(y[[1]]))) )
	) {
		plot <- FALSE
		if (dim>1) {
			warning('Algorithm Successful. Graph Not Easily Plotted Because Values Inhabit Complex Space With Dimension>1.')
		} else {
			all.imaginary <- !any(!Re(y)==0)
			all.real <- !any(!Im(y)==0)
			if (all.real) {
				plot <- TRUE 
			}
			if (all.imaginary) {
				plot(t, Im(y), type='l', las=1, col=4, 
					xlab = 'Domain Of The Solution', 
					ylab = 'Imaginary Part (Simple Affine Linear Spline)', 
					main = 'Approximate Solution (Pure Imaginary)')
			} else {
				plot(y, type='l', las=1, col=4, 
					xlab = 'Real Part Of The Solution Function', 
					ylab = 'Imaginary Part', 
					main = 'Parametric Approximate (Complex) Solution')
				points(as.complex(initial.value))
			}
			abline(v=0, h=0, lty=3)
		}
	}
	if (dim==1) {
		if (calc.dy) {
			results <- data.frame(t=t, y=y, dy=dy)
		} else {
			results <- data.frame(t=t, y=y)
		}
		if (plot) {
			graph <- Cubic.Spline(t, y, dy, method='Hermite', print.all=T, plot=F)$graph
			if (add) {
				lines(graph, col=4)
			} else {
				plot(graph, type='l', las=1, col=4, 
					xlab = 'Domain Of The Solution', 
					ylab = 'Piecewise Hermite Cubic Spline', 
					main = 'Approximate (Spline) Solution To The Diff. Eq.')
			}
			if (points) points(t,y)
			abline(v=0, h=0, lty=3)
		}
	} else {
		if (calc.dy) {
			results <- list(t=t, y=y, dy=dy)
		} else {
			results <- list(t=t, y=y)
		}
		if (plot) {
			if (add) {
				lines(y[[1]], y[[2]], col=3)
				points(initial.value[1], initial.value[2])
			} else {
				ploty1 <- as.numeric(y[[1]])
				ploty2 <- as.numeric(y[[2]])
				plot(ploty1, ploty2, type='l', las=1, col=4, 
					xlab = 'First Coordinate Of Solution Function', 
					ylab = 'Second Coordinate', 
					main = 'Parametric Plot Of The Approximate Solution')
				points(initial.value[1], initial.value[2])
			}
			abline(v=0, h=0, lty=3)
		}
	}
	if (is.null(print)) {if (length(t)<=10) {print <- TRUE} else {print <- FALSE}}
	if (!is.null(store)) assign(store, results, envir = .GlobalEnv)
	if (print) return(results)
#
#	end
#
}



#
##
###
####	Runge-Kutta-Fehlberg
###
##
#

Runge.Kutta.Fehlberg <- function(diff.eq, interval, initial.value, tol=1e-5, min.step.size=1e-4, max.step.size=0.1, plot=TRUE, points=FALSE, print=NULL, store=NULL, add=FALSE, calc.dy=TRUE) {
#
#	warnings
#
Norm <- function(x) sqrt(sum(abs(x)^2))
dim <- length(initial.value)
if (dim>2 & plot) warning('plot=T Will Be Ignored Because dim(initial.value)>2.')
if (!is.function(diff.eq)) stop('f Must Be A Function Of Two Variables, Equal To The Derivative Of Its Second Argument With Respect To Its First Argument, i.e., f(t,y)=dy/dt.')
if (!(is.character(store) | is.null(store))) stop('store Argument Must Be a Character String or NULL.')
if (length(interval)<2) stop('Please Make interval of the Form c(a,b) where a<b.')
if (!interval[2]>interval[1]) stop('Please Make interval of the Form c(a,b) where a<b.')
if (!is.numeric(interval)) stop('Please Make interval of the Form c(a,b) where a<b.')
if (length(interval)>2) warning('interval Has Length Greater Than Two; Places After Two Ignored.')
if (interval[1]+max.step.size > interval[2]) stop('interval Is Too Short (i.e., interval[1]+max.step.size > interval[2]).')
if (!is.numeric(max.step.size) | !max.step.size>0) stop('max.step.size Should Be A Small Positive Real Number.')
if (missing(min.step.size)) {
	min.step.size <- max.step.size/4
} else {
	if (!is.numeric(min.step.size) | !min.step.size>0) stop('min.step.size Should Be A Small Positive Real Number.')
}
stopifnot(min.step.size<max.step.size)
if (!is.logical(plot)) stop('plot Must Be Logical.')
if (!is.logical(points)) stop('points Must Be Logical.')
if (!is.logical(add)) stop('add Must Be Logical.')
if (!(is.logical(print) | is.null(print))) stop('print Should Be Either Logical Or NULL.')
if (tol<=0) stop('tol Must Be Positive (and should be small).')
test <- diff.eq(interval[1], initial.value)
if (!length(test)==dim) stop('The Dimension Of The Codomain Of diff.eq Must Match The Dimension Of initial.vale.')
if (!is.null(dim(test))) stop('We Have !is.null(dim(test)) Where test <- diff.eq(interval[1], initial.value).')
if (!is.numeric(test) & !is.complex(test)) stop('We Have !is.numeric(test) & !is.complex(test) Where test <- diff.eq(interval[1], initial.value).')
#
#	the algorithm
#
	f <- diff.eq
	t <- interval[1]
	inc <- initial.value
	dee <- test
	h <- max.step.size
	i <- go <- 1
	y <- list(inc[1])
	dy <- list(dee[1])
	for (j in 1:dim) {
		y[[j]] <- inc[j]
		dy[[j]] <- dee[j]
		names(y[[j]]) <- paste('y', j, sep='')
		names(dy[[j]]) <- paste('dy', j, sep='')
	}

	while (go==1) {
	flag <- 1
	while (flag==1) {

		K1 <- h*dee
		K2 <- h*f(t[i]+.25*h, inc + .25*K1)
		K3 <- h*f(t[i]+.375*h, inc + .09375*K1 + .28125*K2)
		K4 <- h*f(t[i] + 12*h/13, inc + (1932*K1-7200*K2+7296*K3)/2197)	
		K5 <- h*f(t[i]+h, inc + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104)
		K6 <- h*f(t[i]+0.5*h, inc - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40)

		measure <- Norm( (K1/360 - 128*K3/4275 + 2197*K4/75240 + 0.02*K5 + 2*K6/55)/h )
		
		if (measure<tol & t[i]+h<interval[2]) {
			inc <- inc + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - .2*K5
			for (j in 1:dim) y[[j]][i+1] <- inc[j]
			if (calc.dy) {
				dee <- f(t[i], inc)
				for (j in 1:dim) dy[[j]][i+1] <- dee[j]
			}
			t[i+1] <- t[i]+h
			flag <- 0
			i <- i+1
		} else {
			delta <- 0.84*(tol/measure)^0.25
			if (delta <= 0.1) {h <- 0.1*h} else {if (delta>=4) {h <- 4*h} else {h <- delta*h}}
			if (h>max.step.size) h <- max.step.size
			if (h<min.step.size) h <- min.step.size
			if (t[i]+h > interval[2]) {
				h <- interval[2]-t[i]
				K1 <- h*dee
				K2 <- h*f(t[i]+.25*h, inc + .25*K1)
				K3 <- h*f(t[i]+.375*h, inc + .09375*K1 + .28125*K2)
				K4 <- h*f(t[i] + 12*h/13, inc + (1932*K1-7200*K2+7296*K3)/2197)	
				K5 <- h*f(t[i]+h, inc + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104)
				K6 <- h*f(t[i]+0.5*h, inc - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40)
				inc <- inc + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - .2*K5
				dee <- f(t[i], inc)
				for (j in 1:dim) y[[j]][i+1] <- inc[j]
				if (calc.dy) {
					for (j in 1:dim) dy[[j]][i+1] <- dee[j]
				}
				t[i+1] <- interval[2]
				go <- flag <- 0
				go <- 0
				flag <- 0
				i <- i+1
			}
			if (h==min.step.size) {
				K1 <- h*dee
				K2 <- h*f(t[i]+.25*h, inc + .25*K1)
				K3 <- h*f(t[i]+.375*h, inc + .09375*K1 + .28125*K2)
				K4 <- h*f(t[i] + 12*h/13, inc + (1932*K1-7200*K2+7296*K3)/2197)	
				K5 <- h*f(t[i]+h, inc + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104)
				K6 <- h*f(t[i]+0.5*h, inc - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40)
				inc <- inc + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - .2*K5
				dee <- f(t[i], inc)
				for (j in 1:dim) y[[j]][i+1] <- inc[j]
				if (calc.dy) {
					for (j in 1:dim) dy[[j]][i+1] <- dee[j]
				}
				t[i+1] <- t[i]+h
				flag <- 0
				i <- i+1
			}
		}
	}
	}
	if (dim==1) {
		y <- c(y[[1]])
		if (calc.dy) dy <- c(dy[[1]])
	}
#
#	what to do with the results
#
	if (any(is.complex(y))) {
		plot <- FALSE
		complex.plot <- TRUE
		if (dim>1) {
			warning('Algorithm Successhul. Graph Not Easily Plotted. Values Inhabit Complex Space With Dimension>1.')
		} else {
			all.imaginary <- !any(!Re(y)==0)
			all.real <- !any(!Im(y)==0)
			if (all.real) plot <- TRUE 
			if (all.imaginary) {
				plot(t, Im(y), type='l', las=1, col=4, 
					xlab = 'Domain Of The Solution', 
					ylab = 'Imaginary Part (Simple Affine Linear Spline)', 
					main = 'Approximate Solution (Pure Imaginary)')
			} else {
				plot(y, type='l', las=1, col=4, 
					xlab = 'Real Part Of The Solution Function', 
					ylab = 'Imaginary Part', 
					main = 'Parametric Approximate (Complex) Solution')
				points(as.complex(initial.value))
			}
			abline(v=0, h=0, lty=3)
		}
	}
	if (dim==1) {
		if (calc.dy) {
			results <- data.frame(t=t, y=y, dy=dy)
		} else {
			results <- data.frame(t=t, y=y)
		}
		if (plot) {
			graph <- Cubic.Spline(t, y, dy, method='Hermite', print.all=T, plot=F, table=F)$graph
			if (add) {
				lines(graph, col=4)
			} else {
				plot(graph, type='l', las=1, col=4, 
					xlab = 'Domain Of The Solution', 
					ylab = 'Piecewise Hermite Cubic Spline', 
					main = 'Approximate (Spline) Solution To The Diff. Eq.')
			}
			if (points) points(t,y)
			abline(v=0, h=0, lty=3)
		}
	} else {
		if (calc.dy) {
			results <- list(t=t, y=y, dy=dy)
		} else {
			results <- list(t=t, y=y)
		}
		if (plot) {
			if (add) {
				lines(y[[1]], y[[2]], col=3)
				points(initial.value[1], initial.value[2])
			} else {
				ploty1 <- as.numeric(y[[1]])
				ploty2 <- as.numeric(y[[2]])
				plot(ploty1, ploty2, type='l', las=1, col=4, 
					xlab = 'First Coordinate Of Solution Function', 
					ylab = 'Second Coordinate', 
					main = 'Parametric Plot Of The Approximate Solution')
				points(initial.value[1], initial.value[2])
			}
			abline(v=0, h=0, lty=3)
		}
	}
	if (is.null(print)) {if (length(t)<=10) {print <- TRUE} else {print <- FALSE}}
	if (!is.null(store)) assign(store, results, envir = .GlobalEnv)
	if (print) return(results)
#
#	end
#
}





#
##
###
####		analytically solve a system of two first order ODEs
###
##
#

Solve.Matrix.OIVP <- function(A, initial.time, initial.value) {
#
#	warnings
#
stopifnot(is.matrix(A))
stopifnot(dim(A)[1]==dim(A)[2])

	ei <- eigen(A)
	lam <- ei$values
	v <- ei$vectors
	b <- initial.value
#	if (abs(det(v))<tol) stop('A Does Not Have dim(A)[1]-Many Lin. Indep. Eigenvectors.')
	C <- A
	for (i in 1:dim(A)[1]) {
		for (j in 1:dim(A)[1]) {
			C[i,j] <- v[i,j] * exp( lam[j]*initial.time )
		}
	}
	c <- solve(C,b)
	results <- list(data.frame(scale=c*v[1,], exp=lam))
	for (i in 2:dim(A)[1]) results[[i]] <- data.frame(scale=c*v[i,], exp=lam)
	results

}



Convert.2D <- function(asdf) {

	function (t) {
		summands <- out <- dim <- length(asdf)
		for (i in 1:dim) {
			for (j in 1:dim) {
				summands[j] <- asdf[[i]]$scale[j]*exp( asdf[[i]]$exp[j]*t )
			}
			out[i] <- sum(summands)
		}
		out
#		c( 
#			asdf[[1]]$scale[1]*exp( asdf[[1]]$exp[1]*t ) + asdf[[1]]$scale[2]*exp( asdf[[1]]$exp[2]*t) ,
#			asdf[[2]]$scale[1]*exp( asdf[[2]]$exp[1]*t ) + asdf[[2]]$scale[2]*exp( asdf[[2]]$exp[2]*t)
#		)
	}

}


Solution.Space <- function(f, xlim=c(-5,5), ylim=c(-5,5), buf=0, gridx=NULL, gridy=NULL, add=FALSE, interval=c(0,2)) {

	if (is.null(gridx) | is.null(gridy)) {
		gridspace <- (xlim[2]-xlim[1])/10
		gridspace[2] <- (ylim[2]-ylim[1])/10
		x <- seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/10)
		y <- seq(ylim[1], ylim[2], (ylim[2]-ylim[1])/10)
	} else {
		x <- gridx
		y <- gridy
	}
	xlim <- c(xlim[1]-buf, xlim[2]+buf)
	ylim <- c(ylim[1]-buf, ylim[2]+buf)
	
	if (!add) {
		plot.new()
		par(las=1)
		plot.window(xlim, ylim, las=1)
		axis(1)
		axis(2)
		title(xlab='First Coordinate Of Solution Function(s)')
		title(ylab='Second Coordinate')
		title(main='Scetch Of Two-Dimensional Solution Space')
	} 
	for (i in 1:length(x)) {
		for (j in 1:length(y)) {
			Runge.Kutta(f, interval, c(x[i],y[j]), add=T)
		}
	}
	box()

}


