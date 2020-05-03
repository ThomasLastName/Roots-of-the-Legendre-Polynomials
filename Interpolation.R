
#
##
###	Getting Chebyshev Nodes
##
#


Cheb.Nodes <- function(n, a, b, precise=F, bits=250) {
	nodes <- 1:n
	if (precise) nodes <- mpfr(1:n, bits)
	if (precise) {pie <- Const('pi', bits)} else {pie <- pi}
	for (i in 1:n) {
		nodes[i] <- (a+b)/2 + cos(pie*(2*i-1)/2/n)*(b-a)/2
	}
	rev(nodes)
}

#
##
###	Divided Differences
##
#

Div.Diff <- function(x,y) {

	if (!length(x)==length(y)) {stop('Lengths Of x And y Must Agree.')}
	if (!length(x)==nlevels(factor(x))) {stop('x Contains Repeat Entries. This Does not Work.')}
	if (!is.numeric(x) | !is.numeric(y)) {stop('x And y Must Be Numeric Vectors')}
	if (!is.null(dim(x)) | !is.null(dim(y))) {stop('x And y Must Be Numeric Vectors (not, e.g., matrices).')}
	if (is.list(x) | is.list(y)) stop('x And y Must Be Numeric Vectors')
	if (length(x)<2) {stop('x Too Short.')}

	n <- length(x)
	DivDiff <- diag(n)
	for (j in 1:n) DivDiff[j,1] <- y[j]
	for (j in 2:n) {
		for (i in j:n) {
			DivDiff[i,j] <- (DivDiff[i, j-1] - DivDiff[i-1, j-1]) / (x[i] - x[i-(j-1)])
		}
	}
	DivDiff

}

#
##
###
####	interpolate a value
###
##
#

Interpolate <- function(predictor, response, new, dy=NULL, method='L', precise=F, bits=250) {
	
	if (!length(predictor)==length(response)) {stop('Lengths Of predictor And response Must Agree.')}
	if (!length(predictor)==nlevels(factor(predictor))) {stop('predictor Contains Repeat Entries. This Does not Work.')}
if (!precise) {
	if (is.list(predictor)) stop('predictor And response Must Be Numeric Vectors.')
	if (is.list(response)) stop('predictor And response Must Be Numeric Vectors.')
	if (is.list(new)) stop('new Must Be Numeric Vectors.')
	if (!is.numeric(predictor) | !is.numeric(response)) {stop('predictor And response Must Be Numeric Vectors.')}
	if (!is.null(dim(predictor)) | !is.null(dim(response))) {stop('predictor And response Must Be Numeric Vectors (not, e.g., matrices).')}
	if (!is.numeric(new)) {stop('new (the third argument) Should be a Numeric Vector.')}
	if (!is.null(dim(new))) {stop('new Should be a Numeric Vector.')}	
}
	if (length(predictor)<2) {stop('predictor Too Short.')}
	if (any(!predictor==sort(predictor))) {
		warning('Elements Of response Will Be Used In Asecending Numerical Order, Not In The Order In Which They Are Spceified.')
		response <- response[sort(predictor, index=T)$ix]
		predictor <- sort(predictor)
	}
	stopifnot(method %in% c('B', 'NS', 'CS', 'H', 'L'))
	if (method=='B') warning('Barycentric Representation Is Only Theoretically Valid When Chebyshev Nodes Are Used.')
	
	out <- n <- length(predictor)
	vec <- 1:n
	if (precise) {
		predictor <- mpfr(predictor, bits)
		response <- mpfr(response, bits)
		new <- mpfr(new, bits)
		vec <- mpfr(vec, bits)
		out <- mpfr(out, bits)
		if (!is.null(dy)) dy <- mpfr(dy, bits)
	}
#
#	Lagrange Interpolation
#
	if (method=='L') {
		for (i in 1:length(new)) {
			for (j in 1:n) vec[j] <- prod( (new[i] - predictor[-j]) / (predictor[j] - predictor[-j]) )
			out[i] <- sum(response*vec)
		}
	}
#
#	Hermite Interpolation Based On Derivatives 0 And 1
#
	if (method=='H') {
		if (!precise) {
			if (is.null(dy)) {stop('Hermite Interpolation Requires Data For dy.')}
			if (!is.numeric(dy) | !is.null(dim(dy))) {stop('dy Should Be A Numeric Vector (not, e.g., a matrix).')}
			if (!length(dy)==length(response)) {stop('Lengths Of dy And response Differ.')}
		}

		for (i in 1:length(new)) {
			for (j in 1:n) {
				xx <- predictor[-j]
				vec[j] <- prod( (new[i]-xx)/(predictor[j]-xx) )^2 * ( response[j] - 2*response[j]*(new[i]-predictor[j])*sum(1/(predictor[j]-xx)) + dy[j]*(new[i]-predictor[j]) )
			}
			out[i] <- sum(vec)
		}
	}
#
#	Cubic Spline Interpolation, Natural or Clamped
#
	if (method=='NS' | method=='CS') {
	if (method=='NS') {clamped <- FALSE} else {clamped <- TRUE}
	if (clamped) {
		if (!precise) {
			if (!is.numeric(dy)) {stop('Please Make dy A Vector Of Length 2 Where The First (resp. second) Place Denoted The Derivative Of The Interpolated Function At The Minimum (resp. Maximum) Of The Interval Of Interpolation.')}
			if (!is.null(dim(dy))) {stop('Please Make dy A Vector Of Length 2 Where The First (resp. second) Place Denoted The Derivative Of The Interpolated Function At The Minimum (resp. Maximum) Of The Interval Of Interpolation.')}
			if (!length(dy)==2) {stop('Please Make dy A Vector Of Length 2 Where The First (resp. second) Place Denoted The Derivative Of The Interpolated Function At The Minimum (resp. Maximum) Of The Interval Of Interpolation.')}
		}
		fprimea <- dy[1]
		fprimeb <- dy[2]
	}
	x <- predictor
	y <- response

	thing <- deltac <- deltay <- deltax <- n <- length(x)
	for (i in 1:(n-1)) {deltax[i] <- x[i+1]-x[i]}
	for (i in 1:(n-1)) {deltay[i] <- y[i+1]-y[i]}
	
	if (clamped) {
		A <- list(c(2*deltax[1], deltax[1], rep(0, n-2)))
	} else {
		A <- list(c(1, rep(0, n-1)))
	}
	
		for ( i in 2:(n-1) ) A[[i]] <- c(rep(0, i-2), deltax[i-1], 2*(deltax[i-1] + deltax[i]), deltax[i], rep(0, n-3-(i-2)))
	
	if (clamped) {
		A[[n]] <- c(rep(0, n-2), deltax[n-1], 2*deltax[n-1])
	} else {
		A[[n]] <- c(rep(0, n-1), 1)
	}

	for (i in 1:(n-2)) thing[i] <- deltay[i+1]/deltax[i+1] - deltay[i]/deltax[i]
	if (clamped) 	{
		boundary1 <- deltay[1]/deltax[1] - fprimea
		boundary2 <- fprimeb - deltay[n-1]/deltax[n-1]
	} else {
		boundary1 <- boundary2 <- 0
	}
	thing <- 3*c(boundary1, thing, boundary2)
	
	C <- A[[1]]
	for (i in 2:n) C <- rbind(C, A[[i]])
	soln <- solve(C, thing)
	for (i in 1:(n-1)) {deltac[i] <- soln[i+1]-soln[i]}

	a <- y[1:(n-1)]
	c <- soln[1:(n-1)]
	b <- deltay/deltax - deltax*(deltac + 3*c)/3
	d <- deltac/3/deltax

	for (i in 1:length(new)) {
		z <- new[i]
		j <- max(which(x <= z))
		xj <- x[j]
		out[i] <- a[j] + (z-xj)*( b[j] + (z-xj)*( c[j] + d[j]*(z-xj) ) )
	}
	}
#
#
#
	if (method=='B'){
	w <- (-1)^(1:n-1) * sin( pi*(1:n-0.5)/n )
	C <- (max(x) - min(x))/4
	for (i in 1:length(input)) {
		if (any(input[i]==x)) {
			out[i] <- y[which(input[i]==x)]
		} else {
			out[i] <- sum(w*y/(C*(input[i]-x))) / sum(w/(C*(input[i]-x)))
		}
	}
	}
#
#	print values
#
	if (precise) {
		cbind(new,out)
	} else {
		return( data.frame( New.Data=new, Interpolated.Value=out ) )
	}
#
#	end
#
}





#
##
###
####		Lagrange Interpolation of Data When the Function Is Unknown
###
##
#

Lagrange.Data <- function(x, y, store=NULL, buffer=0, interval=NULL, grain=0.001, method='L', basesgraph=FALSE, plot=TRUE, add=FALSE) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)){
		{stop("store Argument Must Be a Character String or NULL")}
	}
}
if (is.list(x) | is.list(y)) stop('x And y Must Be Numeric Vectors')
if (!nlevels(factor(x))==length(x)) {stop("Nodes Are Not All Distinct.")}
if (!length(x)==length(y)) {stop("Response Data and Predictor Data Are Not Equal in Length.")}
if (length(x)<2) {stop('Too Few Data Points.')}
if (!is.numeric(buffer)) {stop("buffer Argument Must Be Numeric (and non-negative.)")}
if (buffer<0) {stop("buffer Argument Must Be Non-Negative.")}
if (!is.null(dim(x))) {stop("x Must Be A Numeric Vector (not, e.g., a matrix).")}
if (!is.null(dim(y))) {stop("y Must Be A Numeric Vector (not, e.g., a matrix).")}
if (!is.null(interval)) {
	if (!length(interval==2)) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (!interval[2]>interval[1]) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (!is.numeric(interval)) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (length(interval)>2) {warning("interval Has Length More Than Two; Places After Two Ignored.")}
}
if (!is.logical(basesgraph)) {stop('basegraph Must Be Logical.')}
if (!is.logical(plot)) {stop('plot Must Be Logical.')}
if (!is.logical(add)) {stop('add Must Be Logical.')}
if (!( method=='L' | method=='Newton' | method=='Barycentric')) {stop('method May Take One Of Three Values: R, Newton, Barycentric')}
if (method=='Barycentric') warning('Barycentric Representation Is Only Theoretically Valid When Chebyshev Nodes Are Used.')
#
#	the actual algorithm
#
if (is.null(interval)) {
        lo <- min(x) - buffer
        hi <- max(x) + buffer
	} else {
        lo <- interval[1] + buffer
        hi <- interval[2] + buffer
	}
out <- input <- seq(from=lo, to=hi, by=grain)
n <- length(x)
vec <- 1:n
if (method=='L' | basesgraph) {
	plotL <- L <- list(1)
	for (i in 1:n) plotL[[i]] <- plotL[[1]]
	for (i in 1:n) { L[[i]] <- function(z) {prod( (z - x[-i])/(x[i] - x[-i]) )} }
	for (i in 1:n) {for (j in 1:length(input)) { plotL[[i]][j] <- L[[i]](input[j]) } }
	for (i in 1:length(input)) {{for (j in 1:n) vec[j] <- y[j]*plotL[[j]][i]}; out[i] <- sum(vec)}
}
if (method=='Newton') {
	Newton <- diag(Div.Diff(x,y))
	for (i in 1:length(input)) {
		coz <- 1:n
		coz[n] <- Newton[n]
		for (j in (n-1):1) {
			coz[j] <- Newton[j] + coz[j+1]*(input[i] - x[j])
		}
		out[i] <- coz[1]
	}
}
if (method=='Barycentric') {
	w <- (-1)^(1:n-1) * sin( pi*(1:n-0.5)/n )
	C <- (max(x) - min(x))/4
	for (i in 1:length(input)) {
		if (any(input[i]==x)) {
			out[i] <- y[which(input[i]==x)]
		} else {
			out[i] <- sum(w*y/(C*(input[i]-x))) / sum(w/(C*(input[i]-x)))
		}
	}
}
if (any(is.na(out))) warning('NaNs Present Among Interpolated Values')
#
#	Putting Results Into The Global Environment
#
if (!is.null(store)) {
	results <- data.frame(Input=input, Interpolated=out)
	assign(store, results, envir = .GlobalEnv)
}
#
#	graphics
#
if (!add) {
	spreadout <- max(out) - min(out)
	spready <- max(y) - min(y)
	bufx <- (max(x) - min(x))/20
	if (spreadout > 3*spready) {
		upcap <- max(y) + .5*spready
		locap <- min(y) - .5*spready
	}
}
#
#	basesgraph=T
#
if (basesgraph) {
	if (!add) {
		x11()
		if (spreadout > 3*spready) {
			plot(x, y, xlim=c(lo-bufx, hi+bufx), ylim=c(locap, upcap), 
				col=4, las=1, ylab='', main=paste('Weighted Lagrange Basis Polynomials (n=', n,')', sep=''))
		} else {
			plot(x, y, xlim=c(lo-bufx, hi+bufx),
				col=4, las=1, ylab='', main=paste('Weighted Lagrange Basis Polynomials (n=', n,')', sep=''))	
		}		
	} else {
		points(x, y, col=4)
	}
	title(ylab=expression(list(f(x),~~y[1]*L[1](x),~~y[2]*L[2](x),~~ldots,~~y[n]*L[n](x))), mgp=c(2.5,1,0))
	for (i in 1:n) lines(input, y[i]*plotL[[i]], lty=i, col='grey')
}
#
#	plot=T
#
if (plot) {
if (!add) {
	x11()
	if (spreadout > 3*spready) {
		plot(input, out, type='l', col=4, las=1, ylim=c(locap, upcap), 
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Lagrange Interpolation (n=', n, ')', sep='')
			)
	} else {
		plot(input, out, type='l', col=4, las=1,
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Lagrange Interpolation (n=', n, ')', sep='')
			)
	}
	points(x, y)
} else {
	lines(input, out, col=4)
	points(x,y)
}
}
#
#	End
#
}





#
##
###
####		Lagrange Interpolation of Data When the Function Is Known
###
##
#

Lagrange.Auto <- function(f, interval, nodes, method='Barycentric', Chebyshev=TRUE, store=NULL, buffer=0, grain=0.001, basesgraph=FALSE, plot=TRUE, add=FALSE) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)){
		{stop("store Argument Must Be a Character String or NULL")}
	}
}
if (!is.numeric(buffer)) {stop("buffer Argument Must Be Numeric (and non-negative)")}
if (buffer<0) {stop("buffer Argument Must Be Non-Negative")}
if (!nodes==floor(nodes)) {warning("nodes Argument Is Non-Integer; Will Be Rounded Down.")}
if (nodes<2) {stop("Error: Number of Nodes Specified Is Less Than 2")}
if (!is.numeric(nodes)) {stop("nodes Argument Is Not Numeric")}
if (length(nodes)<1) {stop("nodes Argument Has Length Less Than 1")}
if (length(nodes)>1) {warning("length(nodes)>1; Places After 1 Discarded")}
if (length(interval)<2) {stop("interval Argument Must Have Length 2.")}
if (!is.numeric(interval)) {stop("interval Arument Must Be Numeric")}
if (!interval[2]>interval[1]) {stop("Please Make Interval of the Form c(a,b) where a<b.")}
if ( grain > (interval[2]-interval[1])/2 ) {stop("Grain Too Course; Please Set Grain At Most (interval[2]-interval[1])/2.")}
if (!is.logical(basesgraph)) {stop('basegraph Must Be Logical.')}
if (!is.logical(plot)) {stop('plot Must Be Logical.')}
if (!is.logical(add)) {stop('add Must Be Logical.')}
if (!( method=='L' | method=='Newton' | method=='Barycentric')) {stop('method May Take One Of Three Values: L, Newton, Barycentric.')}
if (method=='Barycentric' & Chebyshev==FALSE) warning('Barycentric Representation Is Only Theoretically Valid When Chebyshev Nodes Are Used.')
#
#	algorithm
#
lo <- interval[1]
hi <- interval[2]
buf <- buffer
n <- floor(nodes[1])
if (Chebyshev) {
	x <- Cheb.Nodes(n, lo, hi)
} else {
	x <- seq(from=lo, to=hi, length=n)
}
y <- f(x)
vec <- 1:n
out <- input <- seq(from=lo-buf, to=hi+buf, by=grain)
if (method=='L' | basesgraph) {
	plotL <- L <- list(1)
	for (i in 1:n) plotL[[i]] <- plotL[[1]]
	for (i in 1:n) { L[[i]] <- function(z) {prod( (z - x[-i])/(x[i] - x[-i]) )} }
	for (i in 1:n) {for (j in 1:length(input)) { plotL[[i]][j] <- L[[i]](input[j]) } }
	for (i in 1:length(input)) {{for (j in 1:n) vec[j] <- y[j]*plotL[[j]][i]}; out[i] <- sum(vec)}
}
if (method=='Newton') {
	Newton <- diag(Div.Diff(x,y))
	for (i in 1:length(input)) {
		coz <- 1:n
		coz[n] <- Newton[n]
		for (j in (n-1):1) {
			coz[j] <- Newton[j] + coz[j+1]*(input[i] - x[j])
		}
		out[i] <- coz[1]
	}
}
if (method=='Barycentric') {
	if (!Chebyshev) stop('Cannot Have Both Chebyshev=FALSE and method=Barycentric')
	w <- (-1)^(1:n-1) * sin( pi*(1:n-0.5)/n )
	C <- (max(x) - min(x))/4
	for (i in 1:length(input)) {
		if (any(input[i]==x)) {
			out[i] <- y[which(input[i]==x)]
		} else {
			out[i] <- sum(w*y/(C*(input[i]-x))) / sum(w/(C*(input[i]-x)))
		}
	}
}
if (any(is.na(out))) warning('NaNs Present Among Interpolated Values')
#
#	Putting Results Into The Global Environment
#
if (!is.null(store)) {
	results <- data.frame(Input=input, Interpolated=out)
	assign(store, results, envir = .GlobalEnv)
}
#
#	graphics
#
spreadout <- max(out) - min(out)
spready <- max(y) - min(y)
bufx <- (max(x) - min(x))/20
if (spreadout > 3*spready) {
	upcap <- max(y) + .5*spready
	locap <- min(y) - .5*spready
	}
#
#	basesgraph=T
#
if (basesgraph) {
	if (!add) {
		x11()
		if (spreadout > 3*spready) {
			plot(x, y, xlim=c(lo-bufx, hi+bufx), ylim=c(locap, upcap), 
				col=4, las=1, ylab='', main=paste('Weighted Lagrange Basis Polynomials (n=', n,')', sep=''))
		} else {
			plot(x, y, xlim=c(lo-bufx, hi+bufx),
				col=4, las=1, ylab='', main=paste('Weighted Lagrange Basis Polynomials (n=', n,')', sep=''))	
		}		
	} else {
		points(x, y, col=4)
	}
	title(ylab=expression(list(f(x),~~y[1]*L[1](x),~~y[2]*L[2](x),~~ldots,~~y[n]*L[n](x))), mgp=c(2.5,1,0))
	for (i in 1:n) lines(input, y[i]*plotL[[i]], lty=i, col='grey')
}
#
#	plot=T
#
if (plot) {
if (!add) {
	x11()
	if (spreadout > 3*spready) {
		plot(input, out, type='l', col=4, las=1, ylim=c(locap, upcap), 
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Lagrange Interpolation (n=', n, ')', sep='')
			)
	} else {
		plot(input, out, type='l', col=4, las=1, 
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Lagrange Interpolation (n=', n, ')', sep='')
			)
	}
	points(x, y)
} else {
	lines(input, out, col=4)
	points(x,y)
}
}
#
#	End
#
}





#
##
###
####		Hermite Interpolation Based on Derivatives 0 and 1 When The Function Is Known
###
##
#

Hermite.Auto <- function(f, f1, interval, nodes, Chebyshev=TRUE, store=NULL, points=NULL, buffer=0, grain=0.001, plot=TRUE, add=FALSE) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)) {
		{stop('store Argument Must Be a Character String or NULL.')}
	}
}
if (!is.numeric(buffer)) {stop('buffer Argument Must Be Numeric (and non-negative).')}
if (buffer<0) {stop('buffer Argument Must Be Non-Negative.')}
if (is.null(interval)) {stop('x Is Set to AUTO. So, interval May Not Be Left Unspecified.')}
if (!nodes==floor(nodes)) {warning('nodes Argument Is Non-Integer; Will Be Rounded Down.')}
if (nodes<1) {stop('Error: Number of Nodes Specified Is Less Than 1')}
if (!is.numeric(nodes)) {stop('nodes Argument Is Not Numeric')}
if (length(nodes)<1) {stop('nodes Argument Has Length Less Than 1')}
if (length(nodes)>1) {warning('length(nodes)>1; Places After 1 Discarded')}

	if (length(interval)<2) {stop('Please Make interval of the Form c(a,b) where a<b.')}
	if (!interval[2]>interval[1]) {stop('Please Make interval of the Form c(a,b) where a<b.')}
	if (!is.numeric(interval)) {stop('Please Make interval of the Form c(a,b) where a<b.')}
	if (length(interval)>2) {warning('interval Has Length More Than Two; Places After Two Ignored')}

if (!is.logical(plot)) {stop('plot Must Be Logical.')}
if (!is.logical(add)) {stop('add Must Be Logical.')}
#
#	preulde to the algorithm
#
lo <- interval[1] 
hi <- interval[2] 
buf <- buffer
n <- floor(nodes[1])
if (Chebyshev) {
	x <- Cheb.Nodes(n, lo, hi)
} else {
	x <- seq(from=lo-buf, to=hi+buf, length=n)
}
y <- f(x)
dy <- f1(x)
vec <- 1:n
out <- input <- seq(from=lo-buf, to=hi+buf, by=grain)
#
#	algorithm
#
term <- list(1)
for (i in 1:n) term[[i]] <- term[[1]]
for (i in 1:n) {
	xx <- x[-i]
	for (j in 1:length(input)) {
		term[[i]][j] <- prod( (input[j]-xx)/(x[i]-xx) )^2 * ( f(x[i]) - 2*f(x[i])*(input[j]-x[i])*sum(1/(x[i]-xx)) + f1(x[i])*(input[j]-x[i]) )
	}
}
for (i in 1:length(input)) {
	for (j in 1:n) {
		vec[j] <- term[[j]][i]
	}
	out[i] <- sum(vec[1:n])
}
if (any(is.na(out))) warning('NaNs Present Among Interpolated Values')
#
#	Putting Results Into The Global Environment
#
if (!is.null(store)) {
	results <- data.frame(Input=input, Interpolated=out)
	assign(store, results, envir = .GlobalEnv)
}
#
#	graphics
#
if (plot) {
if (add) {
	lines(input, out, col=4, main=paste('Plot of the Hermite Interpolation (n=', n, ')', sep=''))
	points(x, y)
} else {
	spreadout <- max(out) - min(out)
	spready <- max(y) - min(y)
	if (spreadout > 3*spready) {
		upcap <- max(y) + .5*spready
		locap <- min(y) - .5*spready
	}
	x11()   
	if (spreadout > 3*spready) {
		plot(input, out, type='l', col=4, las=1, ylim=c(locap, upcap),
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Hermite Interpolation (n=', n, ')', sep='')
			)
	} else {
		plot(input, out, type='l', col=4, las=1, 
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Hermite Interpolation (n=', n, ')', sep='')
			)
	}
	points(x, y)
}
}
#
#	End
#
}



#
##
###
####		Hermite Interpolation Based on Derivatives 0 and 1 When The Function Is Generally Unknown
###
##
#

Hermite.Data <- function(x, y, dy, store=NULL, buffer=0, interval=NULL, grain=0.001, plot=TRUE, add=FALSE) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)){
		{stop("store Argument Must Be a Character String or NULL")}
	}
}
if (is.list(x) | is.list(y) | is.list(dy)) stop('x, y, And dy Must Be Numeric Vectors')
if (!is.numeric(x) | !is.numeric(y) | !is.numeric(dy)) {stop("x, y, And dy Must Be Numeric Vectors of Mutually Equal Length")}
if (!nlevels(factor(x))==length(x)) {stop("Nodes Are Not All Distinct.")}
if (!length(x)==length(y)) {stop("x And y Differ In Length.")}
if (!length(dy)==length(y)) {stop("y And dy Differ In Length.")}
if (!length(x)==length(dy)) {stop("x And dy Differ In Length.")}
if (!is.null(dim(x))) {stop("x Must Be A Numeric Vector (not, e.g., a matrix).")}
if (!is.null(dim(y))) {stop("y Must Be A Numeric Vector (not, e.g., a matrix).")}
if (!is.null(dim(dy))) {stop("dy Must Be A Numeric Vector (not, e.g., a matrix).")}
if (!is.numeric(buffer)) {stop("buffer Argument Must Be Numeric (and non-negative)")}
if (buffer<0) {stop("buffer Argument Must Be Non-Negative")}
if (!is.null(interval)) {
	if (length(interval)<2) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (!interval[2]>interval[1]) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (!is.numeric(interval)) {stop("Please Make interval of the Form c(a,b) where a<b.")}
	if (length(interval)>2) {warning("interval Has Length More Than Two; Places After Two Ignored")}
}
if (!is.logical(plot)) {stop('plot Must Be Logical.')}
if (!is.logical(add)) {stop('add Must Be Logical.')}
#
#	prelude to the algorithm
#
if (is.null(interval)) {
        lo <- min(x) - buffer
        hi <- max(x) + buffer
    } else {
        lo <- interval[1] - buffer
        hi <- interval[2] + buffer
}
n <- length(x)
vec <- 1:n
out <- input <- seq(from=lo, to=hi, by=grain)
#
#	algorithm
#
term <- list(1)
for (i in 1:n) term[[i]] <- term[[1]]
for (i in 1:n) {
	xx <- x[-i]
	for (j in 1:length(input)) {
		term[[i]][j] <- prod( (input[j]-xx)/(x[i]-xx) )^2 * ( y[i] - 2*y[i]*(input[j]-x[i])*sum(1/(x[i]-xx)) + dy[i]*(input[j]-x[i]) )
	}
}
for (i in 1:length(input)) {
	for (j in 1:n) {
		vec[j] <- term[[j]][i]
	}
	out[i] <- sum(vec[1:n])
}
if (any(is.na(out))) warning('NaNs Present Among Interpolated Values')
#
#	Putting Results Into The Global Environment
#
if (!is.null(store)) {
	results <- data.frame(Input=input, Interpolated=out)
	assign(store, results, envir = .GlobalEnv)
}
#
#	graphics
#
if (plot) {
if (add) {
	lines(input, out, col=4, main=paste('Plot of the Hermite Interpolation (n=', n, ')', sep=''))
	points(x, y)
} else {
	spreadout <- max(out) - min(out)
	spready <- max(y) - min(y)
	if (spreadout > 3*spready) {
		upcap <- max(y) + .5*spready
		locap <- min(y) - .5*spready
	}
	x11()   
	if (spreadout > 3*spready) {
		plot(input, out, type='l', col=4, las=1, ylim=c(locap, upcap),
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Hermite Interpolation (n=', n, ')', sep='')
			)
	} else {
		plot(input, out, type='l', col=4, las=1, 
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste('Plot of the Hermite Interpolation (n=', n, ')', sep='')
			)
	}
	points(x, y)
}
}
#
#	End
#
}









call.from.frame <- function(x, sto, label=TRUE){
	d2 <- dim(sto)[2]
	stopifnot(is.data.frame(sto))
	stopifnot(d2 > 1)
	if (label) {
		sto[min(which(sto[,1]>=x)),]
	} else {
		if (class(sto[1,2])=='mpfr') {
			k <- mpfr(1:d2, getPrec(sto[1,2]))
		} else {
			k <- 1:d2
		}
		for (i in 1:d2) k[i] <- sto[min(which(sto[,1]>=x)),i]
		k
	}
}



#
##
###
#### 		cubic spline interpolation
###
##
#

Cubic.Spline <- function(x, y, dy=NULL, method='Natural', store=NULL, compute=TRUE, plot=TRUE, points=TRUE, add=FALSE, table=NULL, grain=0.001, print.all=F) {
#
#	warnings
#
	if (!is.numeric(x) | !is.numeric(y)) {stop('x And y Should Be Numeric Vectors.')}
	if (!is.null(dim(x)) | !is.null(dim(y))) {stop('x And y Should Be Numeric Vectors.')}
	if (!nlevels(factor(x))==length(x)) {stop('Elements of x Are Not All Disctinct.')}
	if (is.list(x) | is.list(y)) stop('x And y Must Be Numeric Vectors')	
	if ( !length(y)==length(x) | length(x)<3 ) {stop('x And y Should Be Equal In Length, And Have Length At Least 3.')}
	if (any(!x==sort(x))) {warning('Elements Of x Will Be Used In Asecending Numerical Order, Not In The Order In Which They Are Spceified.')}
		
		x00 <- x
		y00 <- y
		y <- y[sort(x, index=T)$ix]
#		z <- z[sort(x, index=T)$ix]
		x <- sort(x)

	if (!is.logical(plot)) {stop('plot Must Be Logical.')}
	if (!is.logical(print.all)) {stop('print.all Must Be Logical.')}
	if (!is.logical(compute)) {stop('compute Must Be Logical.')}
	if (!is.logical(add)) {stop('add Must Be Logical.')}
	if (!( method=='Natural' | method=='Clamped' | method=='Hermite' )) {stop('method May Take One Of Three Values: Natural, Clamped, Hermite.')}
	if (method=='Clamped') {
		if (!is.numeric(dy)) {stop('Numeric Values Must Be Specified For dy When Method Set To Clamped.')}
		if (!is.null(dim(dy)) | length(dy)<2) {stop('dy Should Be A Numeric Vector of Length 2.')}
		if (length(dy)>2) {warning('dy Has More Than Two Entries; Only The First And Last Entries Will Be Used (for df(a)/dx and df(b)/dx, resp.).')}
	}
	if (method=='Hermite') {
		if (is.null(dy)) {stop('Hermite Interpolation Requires Data For dy.')}
		if (!is.numeric(dy) | !is.null(dim(dy))) {stop('dy Should Be A Numeric Vector (not, e.g., a matrix).')}
		if (!length(dy)==length(x)) {stop('Lengths Of dy And x Differ.')}
		if (is.list(x) | is.list(y) | is.list(dy)) stop('dy Must Be A Numeric Vector')
	}
#	if (!method=='Clamped' & (!is.null(fprimea) | !is.null(fprimeb))) {warning('fprimea And frpimeb Arguments Will Be Ignored Because method Is Not Set To Clamped.')}
	if (method=='Natural' & !is.null(dy)) {warning('dy Argument Will Be Ignored Because method Is Not Set To Hermite.')}
	if (plot & !compute) {
		warning('We Have plot=T Yet compute=F; plot Will Be Ignored, Since There Is Nothing To Plot.')
		plot <- FALSE
	}
	if (compute) {
		if (!is.numeric(grain)) {stop('grain Must Be Numeric')}
		if (!is.null(dim(grain))) {stop('grain Should Be A Numeric Vector Of Length 1 (not, e.g., a matrix.')}
		if (length(grain)<1) {stop('grain Should Be A Small Positive Number.')}
		if (!grain[1]>0) {stop('grain Should Ideally Be A Small Positive Number.')}
		if (length(grain)>1) {warning('length(grain)>1. All But First Place Will Be Ignored.')}
	}
	if (!is.null(table) & !is.logical(table)) {stop('table Must Be Either Logical Or NULL.')}
	if (!is.character(store) & !is.null(store)) {stop("store Argument Must Be a Character String or NULL.")}
#
#	the algorithm
#	
	thing <- deltac <- deltay <- deltax <- n <- length(x)
	for (i in 1:(n-1)) {deltax[i] <- x[i+1]-x[i]}

	if (method=='Natural') {
		Clamped <- FALSE
		Hermite <- FALSE
	}
	if (method=='Clamped') {
		Clamped <- TRUE
		Hermite <- FALSE
		fprimea <- dy[1]
		fprimeb <- dy[length(dy)]
	}
	if (method=='Hermite') {
		Clamped <- FALSE
		Hermite <- TRUE
	}
	if (Hermite) {
#
#	Hermite spline interpolation
#
	a <- b <- c <- d <- 0
	for (i in 1:(n-1)) {
		x0 <- x[i]
		x1 <- x[i+1]
		y0 <- y[i]
		y1 <- y[i+1]
		dy0 <- dy[i]
		dy1 <- dy[i+1]
	sum1 <- c( (x1^2*(y0 - dy0*x0 + (2*x0*y0)/(x0 - x1)))/(x0 - x1)^2, (x1^2*(dy0 - (2*y0)/(x0 - x1)))/(x0 - x1)^2 - (2*x1*(y0 - dy0*x0 + (2*x0*y0)/(x0 - x1)))/(x0 - x1)^2, (y0 - dy0*x0 + (2*x0*y0)/(x0 - x1))/(x0 - x1)^2 - (2*x1*(dy0 - (2*y0)/(x0 - x1)))/(x0 - x1)^2, (dy0 - (2*y0)/(x0 - x1))/(x0 - x1)^2 )
	sum2 <- c(  -(x0^2*(dy1*x1 - y1 + (2*x1*y1)/(x0 - x1)))/(x0 - x1)^2, (2*x0*(dy1*x1 - y1 + (2*x1*y1)/(x0 - x1)))/(x0 - x1)^2 + (x0^2*(dy1 + (2*y1)/(x0 - x1)))/(x0 - x1)^2, - (dy1*x1 - y1 + (2*x1*y1)/(x0 - x1))/(x0 - x1)^2 - (2*x0*(dy1 + (2*y1)/(x0 - x1)))/(x0 - x1)^2, (dy1 + (2*y1)/(x0 - x1))/(x0 - x1)^2 )
		coz <- sum1 + sum2
		a[i] <- coz[1]
		b[i] <- coz[2]
		c[i] <- coz[3]
		d[i] <- coz[4]
	
#		A0 <- y[i]   - x[i]*dy[i]     + 2*x[i]/(x[i]-x[i+1])
#		A1 <- y[i+1] - x[i+1]*dy[i+1] + 2*x[i+1]/(x[i+1]-x[i])
#		I0 <- dy[i] - 2*y[i]/(x[i]-x[i+1])
#		I1 <- dy[i+1] - 2*y[i+1]/(x[i+1]-x[i])

#		a[i] <- A0*( x[i+1]/(x[i+1]-x[i]) )^2 + A1*( x[i]/(x[i+1]-x[i]) )^2
#		b[i] <- I1*( x[i]/(x[i+1]-x[i]) )^2 + I0*( x[i+1]/(x[i+1]-x[i]) )^2 + 2*A1*x[i+1]/(x[i+1]-x[i]) + 2*A0*x[i]/(x[i]-x[i+1])
#		c[i] <- 2*(I0+I1) + (y[i] + y[i+1] - x[i]*dy[i] - x[i+1]*dy[i+1] - 2)/((x[i+1]-x[i])^2)
#		d[i] <- (I0 + I1)/((x[i+1]-x[i])^2)
	}	

	if (compute) {
		lo <- min(x)
		hi <- max(x)
		out <- input <- seq(from=lo, to=hi, by=grain[1])
		for (i in 1:(length(input)-1)) {
			z <- input[i]
			j <- max(which(x <= z))
			xa <- x[j]
			xb <- x[j+1]
			summand1 <- ((z-xb)/(xa-xb))^2 * (y[j] - 2*y[j]*(z-xa)/(xa-xb) + dy[j]*(z-xa))
			summand2 <- ((z-xa)/(xb-xa))^2 * (y[j+1] - 2*y[j+1]*(z-xb)/(xb-xa) + dy[j+1]*(z-xb))
			out[i] <- summand1+summand2
		}
		out[length(input)] <- y[n]
	}

	} else {
#
#	classic spline interploation
#
	for (i in 1:(n-1)) {deltay[i] <- y[i+1]-y[i]}
	
	if (Clamped) {
		A <- list(c(2*deltax[1], deltax[1], rep(0, n-2)))
	} else {
		A <- list(c(1, rep(0, n-1)))
	}
	
		for ( i in 2:(n-1) ) A[[i]] <- c(rep(0, i-2), deltax[i-1], 2*(deltax[i-1] + deltax[i]), deltax[i], rep(0, n-3-(i-2)))
	
	if (Clamped) {
		A[[n]] <- c(rep(0, n-2), deltax[n-1], 2*deltax[n-1])
	} else {
		A[[n]] <- c(rep(0, n-1), 1)
	}

	for (i in 1:(n-2)) thing[i] <- deltay[i+1]/deltax[i+1] - deltay[i]/deltax[i]
	if (Clamped) 	{
		boundary1 <- deltay[1]/deltax[1] - fprimea
		boundary2 <- fprimeb - deltay[n-1]/deltax[n-1]
	} else {
		boundary1 <- boundary2 <- 0
	}
	thing <- 3*c(boundary1, thing, boundary2)
	
	C <- A[[1]]
	for (i in 2:n) C <- rbind(C, A[[i]])
	soln <- solve(C, thing)
	for (i in 1:(n-1)) {deltac[i] <- soln[i+1]-soln[i]}

	a <- y[1:(n-1)]
	c <- soln[1:(n-1)]
	b <- deltay/deltax - deltax*(deltac + 3*c)/3
	d <- deltac/3/deltax
	
	if (compute) {
		lo <- min(x)
		hi <- max(x)
		out <- input <- seq(from=lo, to=hi, by=grain[1])
		for (i in 1:(length(input)-1)) {
			z <- input[i]
			j <- max(which(x <= z))
			xj <- x[j]
			out[i] <- a[j] + (z-xj)*( b[j] + (z-xj)*( c[j] + d[j]*(z-xj) ) )
		}
		out[length(input)] <- y[n]
		if (any(is.na(out))) warning('NaNs Present Among Interpolated Values')
	} else {
		out <- input <- NA
	}
	
	}
#
#	plot = T
#
if (plot) {
	if (!add) {
		plot(input, out, las=1, col=4, type='l',
			xlab = 'Region of Interpolation',
			ylab = 'Interpolated Value',
			main = paste(
					method,
					' Spline Interpolant Based On ',
					length(x),
					' Points',
				sep=''
				)
			)
		} else {
			lines(input, out, col=4)
		}
	if (is.null(points)) points <- TRUE
	if (points) points(x00,y00)
}
#
#	what to do with the results
#
	results <- list(
		meta=c('Letting x_0,...,x_n denote the location of the given nodes,', 
			paste('the', method, ' cubic spline interpolant is of the form'), 
			'S(x)= a_{0,i} + a_{1,i}*(x-x_i) + a_{2,i}*(x-x_i)^2 + a_{3,i}*(x-x_i)^3))',
			'for i= 0,1,...,(length(x)-2), with coefficients depending on the interval',
			'in which x lies, as specified by the following table.'
			),
		coefficients=data.frame( 
			i=0:(n-2),
			x_i=x[1:(n-1)],
			When_x_In=paste('[',x[1:(n-1)], ' ', x[1:(n-1)] + deltax, ']', sep=''),
			Affine.Coeff=a, 
			Linear.Coeff=b, 
			Quadratic.Coeff=c, 
			Cubic.Coeff=d
			),
		concise=as.matrix(data.frame(
			Affine=a,
			Linear=b,
			Quadratic=c,
			Cubic=d
			)),
		graph=data.frame(
			Input=input,
			Interpolated=out
			),
		nodes=data.frame(
			Location_of_Nodes=x00,
			Value_at_Nodes=y00
			)
		)
	
	if (!is.null(store)) assign(store, results, envir=.GlobalEnv)
#
#
#
	if (is.null(table)) {
		if (n>10) {
			table <- FALSE
		} else {
			table <- TRUE
		}
	}
	if (table) print(results$coefficients)
	if (print.all) results
#
#	end
#
}


