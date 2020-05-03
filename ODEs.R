




#
##
###	implement an adaptive step-size scheme in an explicit one step method
##
#

Adapt.One.Step.Method <- function(Increment.Function, Initial.Step.Size, Starting.Point, Starting.Value, tol, K, m, Length.Of.Domain) {

#	m is the order of consistency of the one step method with increment function Increment.Function

	h <- Initial.Step.Size
	t <- Starting.Point
	y <- Starting.Value
	I <- Increment.Function		# function of the form f(t[i-1], t[i], diff.eq, y[i-1])
	T <- Length.Of.Domain
	go <- 1
	
	while(go==1) {
		go <- 0
		H <- 2*h
		one.big <- I(t, t+H, f, y)
		two.small <- I(t, t+h, f, y)
		two.small <- I(t+h, t+H, f, two.small)
		C <- abs( (two.small-one.big) / (H^(m+1)*(1-2^(-m))) )
		h <- (tol/(K*T*C))^(1/m)
		if (h<H/4) go <- 1
	}
	
	(2^m * two.small - one.big) / (2^m-1)

}





#
##
###	the increment function in Heun's method
##
#

Increment.Heun <- function(from, to, f, y.at.from) {

	thing <- f(from, y.at.from)
	y.at.from + 0.5*(to-from)*( thing + f(to, y.at.from+(to-from)*thing) )

}

Step <- function(h, t, y) Adapt.One.Step.Method <- function(Increment.Heun, h, t, y, 1e-5, 10, 2, 2.9)









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

RK4.Test <- function(diff.eq, interval, initial.value, step.size) {

	start <- Sys.time()
	f <- diff.eq
#	t <- seq(interval[1], interval[2], step.size)
	t <- leng <- floor( (interval[2]-interval[1])/step.size )
	for (x in 1:leng) t[x] <- interval[1] + x*(interval[2]-interval[1])/(leng-1)
	y <- initial.value
	for (i in 2:length(t)) y[i] <- Increment.RK4(t[i-1], t[i], f, y[i-1])
	plot(t, y, las=1, col=4, type='l', xlab='Time Domain', ylab='Approximate Value of the Solution Function', main='Results of RK4')
	print(paste( 'Time Required for Computation: ', Sys.time()-start ))
}

RK4.Test(f, c(-3,3), 1/1901, 0.000001)


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

plot.2D <- function(g, interval, grain=0.005, add=FALSE) {

	x <- y <- input <- seq(interval[1], interval[2], grain)
	for (i in 1:length(input)) {
		x[i] <- g(input[i])[1]
		y[i] <- g(input[i])[2]
	}
	if (add) {
		lines(x, y, lty=3)
	} else {
		plot(x, y, type='l', las=1)
	}
	lines(v=0, h=0, lty=2)
}

Norm <- function(x) sqrt(sum(abs(x)^2))

Runge.Kutta.Fehlberg(f, c(0,2*pi), c(1,1))


Mat.Test <- function(A, interval, initial.value, custom.h=exp(-seq(2, 12, 0.1)), method='RK') {

	f <- function(t,y) c( A %*% y )
	asdf <- Solve.Matrix.OIVP(A, interval[1], initial.value)
	ytrue <- Convert.2D(asdf)
	Test(f, ytrue, interval, initial.value, custom.h, method)

}

Test <- function(f, ytrue, interval, initial.value, custom.h=exp(-seq(2, 12, 0.1)), method='RK') {
	
	h <- custom.h
	dim <- length(initial.value)
	control <- results <- list(0)
	count <- compgood <- compbad <- co <- comparison <- err <- error <- numeric(0)
	for (i in 1:length(h)) {
		if (method=='RK') {
			results[[i]] <- Runge.Kutta(f, interval, initial.value, h[i], plot=F, print=T, calc.dy=FALSE)
		} else {
			results[[i]] <- Euler.1(f, interval, initial.value, h[i], dim)
		}
		control[[i]] <- list(0)
		for (j in 1:dim) {
			control[[i]][[j]] <- 0
			for (k in 1:length(results[[i]]$t)) {
				control[[i]][[j]][k] <- ytrue( results[[i]]$t[k] )[j] 
			}
		}
		for (k in 1:length(results[[i]]$t)) {
			errrr <- numeric(dim)
			for (j in 1:dim) errrr[j] <- control[[i]][[j]][k] - results[[i]]$y[[j]][k]
			err[k] <- Norm( errrr )
		}
		error[i] <- max(err)
	}

	complete <- data.frame(h=h, Nh=error)
	one.step <- (log(error[2:101]) - log(error[1:100])) / (log(h[2:101]) - log(h[1:100]))
	two.step <- ((one.step[2:100]) - (one.step[1:99])) / (log(h[3:101]) - log(h[1:99]))
	two.stepp <- two.step
	two.stepp[1] <- 2
	for (i in min(which(abs(two.stepp)>1)):99) count[i] <- i-max(which(abs(two.stepp[min(which(abs(two.stepp)>1)):i])>1))
#	collapse <- which(count==floor(max(count)/2)) + min(which(abs(two.step[which(count==floor(max(count)/2)):99])>1))
	collapse <- 1+max(which(count>2))


	hgood <- h[1:collapse]
	hbad <- h[collapse:101]
	errgood <- error[1:collapse]
	errbad <- error[collapse:101]

if (FALSE) {
	eta <- seq(1, 8, 0.01)
	for (i in 1:length(eta)) {
		et <- eta[i]
		co[i] <- lsfit(h^et, error)$coefficients[2]
		comparison[i] <- sum( ( error-co[i]*h^et )^2 )
		M <- co
	}
}

#	if (length(hgood)>0) {
#		for (i in 1:length(eta)) {
#			et <- eta[i]
#			co[i] <- lsfit(hgood^et, errgood)$coefficients[2]
#			compgood[i] <- sum( ( errgood-co[i]*hgood^et )^2 )
#		}
#	} else {
#		compgood <- NA
#	}

#	if (length(hbad)>0) {
#		for (i in 1:length(eta)) {
#			et <- eta[i]
#			co[i] <- lsfit(hbad^et, errbad)$coefficients[2]
#			compbad[i] <- sum( ( errbad-co[i]*hbad^et )^2 )
#		}
#	} else  {
#		compbad <- NA
#	}
	
#	plot(two.step, las=1, main=expression(Discrete~Approximation~Of~d^2~log(N(h))/d~log(h)^2))
#	abline(v=collapse)
if (FALSE) {	
	summary <- c( eigen(A)$values,  
			as.numeric( collapse ),
			as.numeric( h[collapse] ),
			as.numeric( interval[1] ),
			initial.value,
			as.numeric( eta[which(comparison==min(comparison))] ),
			as.numeric( M[which(comparison==min(comparison))] ),
#			eta[which(compgood==min(compgood))],
#			eta[which(compbad==min(compbad))],
			as.numeric( mean( one.step ) ), 
			as.numeric( mean( one.step[1:collapse] ) ),
			as.numeric( mean( one.step[collapse:100] ) ),
			as.numeric( min(comparison) )
		)
	names(summary) <- c(
			'Eigen1', 
			'Eigen2',
			'Collapse',
			'Size Of Collapse', 
			't0',
			'y(t0)_1',
			'y(t0)_2',
			'eta hat 1',
			'Fitted',
#			'eta hat 1 pre-collapse',
#			'eta hat 1 post-collapse',
			'eta hat 2', 
			'eta hat 2 pre-collapse',
			'eta hat 2 post-collapse',
			'Sum of Squares'
		)
}
	
	summary <- list(
			matrix = A,
			eigen = eigen(A, only.values=TRUE)$values,
			Collapse = as.numeric( collapse ),
			Size.Of.Collapse = as.numeric( h[collapse] )
		)
	list(
		summary = summary,
		complete = data.frame(h=h, Nh=error)
	)

	#	eig1, eig2, eta hat, alt eta hat, summed squared, a11, a12, a22, a22

}


Discrete.Derivative <- function(x,y) {
	stopifnot(length(x)==length(y))
	stopifnot(length(x)>2)
	M <- 2:length(x)
	m <- 1:(length(x)-1)
	(y[M]-y[m]) / (x[M]-x[m])
}


DSD <- function(x,y) {
	dy <- Discrete.Derivative(x,y)
	( dy[2:(length(x)-1)]-dy[1:(length(x)-2)] ) / ( x[3:length(x)]- x[1:(length(x)-2)])
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



Re_Test <- function(stuff, new.h) {

	interval <- c(0,0.5)
	A <- stuff$summary$matrix
	initial.value <- rep(1, dim(A)[1])
	h <- new.h
	f <- function(t,y) c( A %*% y )
	asdf <- Solve.Matrix.OIVP(A, interval[1], initial.value)
	ytrue <- Convert.2D(asdf)

	control <- results <- list(0)
	count <- compgood <- compbad <- co <- comparison <- err <- error <- numeric(0)
	for (i in 1:length(h)) {
		results[[i]] <- Runge.Kutta(f, interval, initial.value, h[i], plot=F, print=T)
		control[[i]] <- list(0)
		for (j in 1:2) {
			control[[i]][[j]] <- 0
			for (k in 1:length(results[[i]]$t)) {
				control[[i]][[j]][k] <- ytrue( results[[i]]$t[k] )[j] 
			}
		}
		for (k in 1:length(results[[i]]$t)) err[k] <- Norm( c(control[[i]][[1]][k] - results[[i]]$y[[1]][k], control[[i]][[2]][k] - results[[i]]$y[[2]][k]) )
		error[i] <- max(err)
	}
	h <- c(stuff$complete$h, h)
	error <- c(stuff$complete$Nh, error)

	complete <- data.frame(h=h, Nh=error)
	one.step <- (log(error[2:101]) - log(error[1:100])) / (log(h[2:101]) - log(h[1:100]))
	two.step <- ((one.step[2:100]) - (one.step[1:99])) / (log(h[3:101]) - log(h[1:99]))
	two.stepp <- two.step
	two.stepp[1] <- 2
	for (i in min(which(abs(two.stepp)>1)):99) count[i] <- i-max(which(abs(two.stepp[min(which(abs(two.stepp)>1)):i])>1))
#	collapse <- which(count==floor(max(count)/2)) + min(which(abs(two.step[which(count==floor(max(count)/2)):99])>1))
	collapse <- 1+max(which(count>2))
	
	summary <- list(
			matrix = A,
			eigen = eigen(A, only.values=TRUE)$values,
			Collapse = as.numeric( collapse ),
			Size.Of.Collapse = as.numeric( h[collapse] )
		)
	list(
		summary = summary,
		complete = data.frame(h=h, Nh=error)
	)

}








Robust.Test <- function(A, interval, initial.value, custom.h=exp(-seq(2, 12, 0.1))) {
	
	h <- custom.h
	f <- function(t,y) c( A %*% y )
	asdf <- Solve.Matrix.OIVP(A, interval[1], initial.value)
	ytrue <- Convert.2D(asdf)

#	control <- results <- list(0)
	error <- rep(0, length(h))
	count <- compgood <- compbad <- co <- comparison <- err <- numeric(0)

	for (i in 1:length(h)) {
		leng <- floor( (interval[2]-interval[1])/h[i] )
		t <- c( interval[1], interval[1]+h )
		inc <- initial.value
		for (l in 2:leng) {						# l not 1
			inc <- Increment.RK4(t[1], t[2], f, inc)		# by contrast, this is 1
			incinc <- Norm( ytrue(t[2]) - inc )
			if (incinc>error[i]) error[i] <- incinc
			t <- c( t[2], t[2]+h )
		}
	}

	complete <- data.frame(h=h, Nh=error)
	one.step <- (log(error[2:101]) - log(error[1:100])) / (log(h[2:101]) - log(h[1:100]))
	two.step <- ((one.step[2:100]) - (one.step[1:99])) / (log(h[3:101]) - log(h[1:99]))
	two.stepp <- two.step
	two.stepp[1] <- 2
	for (i in min(which(abs(two.stepp)>1)):99) count[i] <- i-max(which(abs(two.stepp[min(which(abs(two.stepp)>1)):i])>1))
#	collapse <- which(count==floor(max(count)/2)) + min(which(abs(two.step[which(count==floor(max(count)/2)):99])>1))
	collapse <- 1+max(which(count>2))


	hgood <- h[1:collapse]
	hbad <- h[collapse:101]
	errgood <- error[1:collapse]
	errbad <- error[collapse:101]

if (FALSE) {
	eta <- seq(1, 8, 0.01)
	for (i in 1:length(eta)) {
		et <- eta[i]
		co[i] <- lsfit(h^et, error)$coefficients[2]
		comparison[i] <- sum( ( error-co[i]*h^et )^2 )
		M <- co
	}
}

#	if (length(hgood)>0) {
#		for (i in 1:length(eta)) {
#			et <- eta[i]
#			co[i] <- lsfit(hgood^et, errgood)$coefficients[2]
#			compgood[i] <- sum( ( errgood-co[i]*hgood^et )^2 )
#		}
#	} else {
#		compgood <- NA
#	}

#	if (length(hbad)>0) {
#		for (i in 1:length(eta)) {
#			et <- eta[i]
#			co[i] <- lsfit(hbad^et, errbad)$coefficients[2]
#			compbad[i] <- sum( ( errbad-co[i]*hbad^et )^2 )
#		}
#	} else  {
#		compbad <- NA
#	}
	
#	plot(two.step, las=1, main=expression(Discrete~Approximation~Of~d^2~log(N(h))/d~log(h)^2))
#	abline(v=collapse)
if (FALSE) {	
	summary <- c( eigen(A)$values,  
			as.numeric( collapse ),
			as.numeric( h[collapse] ),
			as.numeric( interval[1] ),
			initial.value,
			as.numeric( eta[which(comparison==min(comparison))] ),
			as.numeric( M[which(comparison==min(comparison))] ),
#			eta[which(compgood==min(compgood))],
#			eta[which(compbad==min(compbad))],
			as.numeric( mean( one.step ) ), 
			as.numeric( mean( one.step[1:collapse] ) ),
			as.numeric( mean( one.step[collapse:100] ) ),
			as.numeric( min(comparison) )
		)
	names(summary) <- c(
			'Eigen1', 
			'Eigen2',
			'Collapse',
			'Size Of Collapse', 
			't0',
			'y(t0)_1',
			'y(t0)_2',
			'eta hat 1',
			'Fitted',
#			'eta hat 1 pre-collapse',
#			'eta hat 1 post-collapse',
			'eta hat 2', 
			'eta hat 2 pre-collapse',
			'eta hat 2 post-collapse',
			'Sum of Squares'
		)
}
	
	summary <- list(
			matrix = A,
			eigen = eigen(A, only.values=TRUE)$values,
			Collapse = as.numeric( collapse ),
			Size.Of.Collapse = as.numeric( h[collapse] )
		)
	list(
		summary = summary,
		complete = data.frame(h=h, Nh=error)
	)

	#	eig1, eig2, eta hat, alt eta hat, summed squared, a11, a12, a22, a22

}

Robust.Test(A, c(0,0.5), c(1,1), custom.h=1e-3)




Euler.1 <- function(diff.eq, interval, initial.value, step.size, dim) {

	lo <- interval[1]
	hi <- interval[2]
	h <- step.size
	t <- seq(lo, hi, h)
	inc <- initial.value
	y <- list(inc[1])
	for (i in 1:dim) {
		y[[i]] <- inc[i]
		names(y[[i]]) <- paste('y', i, sep='')
	}
	for (i in 2:length(t)) {
		inc <- inc + (t[i]-t[i-1])*f(t[i-1], inc)
		for (j in 1:dim) y[[j]][i] <- inc[j]
	}
	
	list(t=t, y=y)

}





f <- function(t,y) c(
	cos(y[1]) - cos(y[2]),
	sin(y[1]) - sin(y[2]) + t
)

if (TRUE) {
	T <- 10
	buf <- 1
	xlim <- c(-20,20)
	ylim <- c(-10,10)
	x <- seq(xlim[1]-buf, xlim[2]+buf, len=40)
	y <- seq(ylim[1]-buf, ylim[2]+buf, len=20)
	if (TRUE) {
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
			Runge.Kutta(f, interval=c(0,T), c(x[i],y[j]), step=1e-2, plot=F, store='sto')
			lines(sto$y[[1]], sto$y[[2]], col=3)
#			points(sto$y[[1]][1], sto$y[[2]][1], cex=0.3)
		}
	}
	box()
#	for (i in 0:20) {abline(2*pi*i,1, lty=3); abline(-2*pi*i,1, lty=3)}
}




