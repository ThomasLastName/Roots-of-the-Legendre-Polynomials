

#
##
###	the good 'ol bisection method for root finding
##
#

Bisection.Method <- function(f, interval, times=100, store=NULL, plot=TRUE, lines=TRUE, print=NULL, check=TRUE, add=FALSE, precise=F, bits=250) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)) {stop("store Argument Must Be a Character String or NULL")}
}
if (!is.function(f)) {stop('f Must Be A Function Of x')}
if (length(interval)<2) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (!is.numeric(interval)) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (!interval[2]>interval[1]) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (length(interval)>2) {warning("interval Has Length More Than Two; Places After Two Ignored")}
if (!is.numeric(times)) {stop("times Must Be A Positive Number")}
if (!is.null(dim(times))) {stop("times Must Be A Positive Number (not, e.g., a matrix)")}
if (times<1) {stop("times Must Be A Positive Number")}
if (length(times)>1) {warning("times Has Length Greater Than 1, Only The First Element Will Be Used")}
if (!floor(times[1])==times[1]) {warning("times Is Non-Integer, Will Be Rounded Down")}

if (!is.numeric(bits)) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (!is.null(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (is.list(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (length(bits)>1) warning('Only The First Entry Of bits Will Be Used.')
if (bits<1) stop('bits Must Be At Least 1 (and should probably be at least 53).')
if (!is.logical(precise)) stop('precise Must Be Logical.')	
if (precise & !'Rmpfr'%in%(.packages())) stop('In Order To Utilize The Setting precise=T, The Package Rmpfr Is Required.')
#
#	the actual algorithm
#
if (check) {
	trial <- seq(from=interval[1], to=interval[2], by=0.001)
	trial <- sign(f(trial))
	if (length(unique(trial))==1) {
		warning('Seems Like f Might Not Have A Root In The Specified Interval.')
	}
}
if (precise) {
	a <- b <- c <- mpfr(1:(times+1), bits)           # Create three numeric vectors of length 21=(initial guess)+(20 iterations)
	interval <- mpfr(interval, bits)
} else {
	a <- b <- c <- 1:(times+1)                       # Create three numeric vectors of length 21=(initial guess)+(20 iterations)
}
a[1] <- interval[1]                                    # Designate initial value
b[1] <- interval[2]                                    # Designate initial value
for (i in 1:(times+1)){
        c[i] <- a[i] + (b[i]-a[i])/2;
        if ( sign(f(c[i]))*sign(f(a[i]))>0 ) {a[i+1]<-c[i]; b[i+1]<-b[i]} 
        else {a[i+1]<-a[i]; b[i+1]<-c[i]}
    }
if ( any(f(c)==0) ) for(i in min(which(f(c)==0)):(times+1)) {                 # Correct in case of accidental success
        a[i] <- a[min(which(f(c)==0))]; 
        b[i] <- b[min(which(f(c)==0))]; 
        c[i] <- c[min(which(f(c)==0))]
    }
a <- a[1:(times+1)]; b <- b[1:(times+1)];                                     # Cut off vectors a and b after 21 places
#
#	what to do with the results
#
if (precise) {
	results <- c
} else {
	results <- data.frame(Iteration=0:times, Lower=a, Upper=b, Estimated_Root=c)
}
if (!is.null(store)) {                                                        # Store results in a data frame inthe gloabl environment
	assign(store, results, envir = .GlobalEnv)
}
if (is.null(print)) {                                                         # Show results in a table
	if (times<16) {
		print(results)
	} else {   
		print(tail(results))
	}
}
if (add) {
if (plot) {                                                                   # Show results in a graph
	points(0:times, c, las=1, pch=4,
		xlab='Iterations', 
		ylab='Estimated Root', 
		main='Results of The Bisection Method')
	if (lines) {for (i in 1:(times+1)) segments(i-1, as.numeric(a[i]), i-1, as.numeric(b[i]), lty=2)}
	}
} else{
if (plot) {
	x11()                                                                   # Show results in a graph
	plot(0:times, c, las=1, 
		xlab='Iterations', 
		ylab='Estimated Root', 
		main='Results of The Bisection Method')
	if (lines) {for (i in 1:(times+1)) segments(i-1, as.numeric(a[i]), i-1, as.numeric(b[i]), lty=2)}
	}
}
if (is.logical(print)){
	if (print) results
}
#
#	end
#
}



#
##
###	regula falsi, otherwise known as the method of false position
##
#

Regula.Falsi <- function(f, interval, times=100, store=NULL, plot=TRUE, lines=TRUE, print=NULL, check=TRUE, add=FALSE, tol=NULL, precise=F, bits=250) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)) {stop("store Argument Must Be a Character String or NULL")}
}
if (!is.function(f)) {stop('f Must Be A Function Of x')}
if (!is.null(tol)) {
	if (!is.numeric(tol)) {stop('tol Must Be Positive Real (and should be very small)')}
	if (tol<=0) {stop('tol Must Be Positive Real (and should be very small)')}
} else {
	if (precise) {tol <- 1e-10} else {tol <- 1e-30}
}
if (length(interval)<2) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (!is.numeric(interval)) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (!interval[2]>interval[1]) {stop("Please Make interval of the Form c(a,b) where a<b.")}
if (length(interval)>2) {warning("interval Has Length More Than Two; Places After Two Ignored")}
if (!is.numeric(times)) {stop("times Must Be A Positive Number")}
if (!is.null(dim(times))) {stop("times Must Be A Positive Number (not, e.g., a matrix)")}
if (times<1) {stop("times Must Be A Positive Number")}
if (length(times)>1) {warning("times Has Length Greater Than 1, Only The First Element Will Be Used")}
if (!floor(times[1])==times[1]) {warning("times Is Non-Integer, Will Be Rounded Down")}

if (!is.numeric(bits)) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (!is.null(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (is.list(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (length(bits)>1) warning('Only The First Entry Of bits Will Be Used.')
if (bits<1) stop('bits Must Be At Least 1 (and should probably be at least 53).')
if (!is.logical(precise)) stop('precise Must Be Logical.')	
if (precise & !'Rmpfr'%in%(.packages())) stop('In Order To Utilize The Setting precise=T, The Package Rmpfr Is Required.')
#
#	the actual algorithm
#
if (check) {
	trial <- seq(from=interval[1], to=interval[2], by=0.001)
	trial <- sign(f(trial))
	if (length(unique(trial))==1) {
		warning('Seems Like f Might Not Have A Root In The Specified Interval.')
	}
}
if (precise) {
	a <- b <- c <- mpfr(1:(times+1), bits)           # Create three numeric vectors of length 21=(initial guess)+(20 iterations)
	interval <- mpfr(interval, bits)
} else {
	a <- b <- c <- 1:(times+1)                       # Create three numeric vectors of length 21=(initial guess)+(20 iterations)
}
a[1] <- interval[1]                                    # Designate initial value
b[1] <- interval[2]                                    # Designate initial value
for (i in 1:(times+1)){                                # The main event
        if ( any(abs(f(a[1:i])-f(b[1:i]))<tol) ) {c[i] <- 666} else { 
        c[i] <- (b[i]*f(a[i]) - a[i]*f(b[i]))/( f(a[i])-f(b[i]) )}
        if ( sign(f(c[i]))*sign(f(a[i]))>0 ) {a[i+1]<-c[i]; b[i+1]<-b[i]} 
        else {a[i+1]<-a[i]; b[i+1]<-c[i]}
    }
if ( any(f(c)==0) ) for(i in min(which(abs(f(c))<tol)):(times+1)) {           # Correct in case of accidental success
        a[i] <- a[min(which(abs(f(c))<tol))]; 
        b[i] <- b[min(which(abs(f(c))<tol))]; 
        c[i] <- c[min(which(abs(f(c))<tol))]
    }
a <- a[1:(times+1)]; b <- b[1:(times+1)];  
#
#	what to do with the results
#
if (precise) {
	results <- c
} else {
	results <- data.frame(Iteration=0:times, Lower=a, Upper=b, Estimated_Root=c)
}
if (!is.null(store)) {                                                        # Store results in a data frame inthe gloabl environment
	assign(store, results, envir = .GlobalEnv)
}
if (is.null(print) & !any(c==666)) {                                                         # Show results in a table
	if (times<16) {
		print(results)
	} else {   
		print(tail(results))
	}
}
if (add) {
if (plot & !any(c==666)) {                                                                   # Show results in a graph
	points(0:times, c, las=1, pch=4,
		xlab='Iterations', 
		ylab='Estimated Root', 
		main='Results of The Regula Falsi')
	if (lines) {for (i in 1:(times+1)) segments(i-1, a[i], i-1, b[i], lty=2)}
	}
} else {
if (plot & !any(c==666)) {
	x11()                                                                   # Show results in a graph
	plot(0:times, c, las=1, 
		xlab='Iterations', 
		ylab='Estimated Root', 
		main='Results of The Regula Falsi')
	if (lines) {for (i in 1:(times+1)) segments(i-1, a[i], i-1, b[i], lty=2)}
	}
}
if (any(c==666)) {                                                            # a half-assed attempt to put safety on a foot gun
	print(results[min(which(c==666)),])
	print(paste('We had f(a_n)=f(b_n) at some point, indicated by 666'))
	}
if (is.logical(print)){
	if (print) results
}
#
#	end
#

}




#
##
###	secant method
##
#

Secant.Method <- function(f, guess1, guess2=NULL, times=100, store=NULL, plot=TRUE, print=NULL, check=TRUE, add=FALSE, tol=NULL, precise=F, bits=250) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)) {stop("store Argument Must Be a Character String or NULL")}
}
if (!is.function(f)) {stop('f Must Be A Function Of x')}
if (!is.null(tol)) {
	if (!is.numeric(tol)) {stop('tol Must Be Positive Real (and should be very small)')}
	if (tol<=0) {stop('tol Must Be Positive Real (and should be very small)')}
} else {
	if (precise) {tol <- 1e-10} else {tol <- 1e-30}
}
if (!is.numeric(guess1)) {stop("Please Make First Initial Guess a Real Number")}
if (!is.null(dim(guess1))) {stop("Please Make First Initial Guess a Real Number (not, e.g., a matrix)")}
if (length(guess1)>1) {warning("guess1 Has Length Greater Than 1, Only The First Element Will Be Used")}

if (is.null(guess2)) {return('Secant Method Requires Two Initial Guesses')}
if (!is.numeric(guess2)) {stop("Please Make Second Initial Guess a Real Number")}
if (!is.null(dim(guess2))) {stop("Please Make Second Initial Guess a Real Number (not, e.g., a matrix)")}
if (length(guess2)>1) {warning("guess2 Has Length Greater Than 1, Only The First Element Will Be Used")}

if (guess1==guess2) {stop('Initial Guesses Must Be Distinct')}
if (!is.numeric(times)) {stop("times Must Be A Positive Number")}
if (!is.null(dim(times))) {stop("times Must Be A Positive Number (not, e.g., a matrix)")}
if (times<1) {stop("times Must Be A Positive Number")}
if (length(times)>1) {warning("times Has Length Greater Than 1, Only The First Element Will Be Used")}
if (!floor(times[1])==times[1]) {warning("times Is Non-Integer, Will Be Rounded Down")}

if (!is.numeric(bits)) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (!is.null(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (is.list(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (length(bits)>1) warning('Only The First Entry Of bits Will Be Used.')
if (bits<1) stop('bits Must Be At Least 1 (and should probably be at least 53).')
if (!is.logical(precise)) stop('precise Must Be Logical.')	
if (precise & !'Rmpfr'%in%(.packages())) stop('In Order To Utilize The Setting precise=T, The Package Rmpfr Is Required.')
#
#	the actual algorithm
#
if (check) {
	ey <- min(guess1, guess2)
	bee <- max(guess1, guess2)
	ey0 <- ey - min(0.2, (bee-ey)/2)
	bee0 <- bee + min(0.2, (bee-ey)/2)
	trial <- seq(ey0, bee0, by=0.001)
	trial <- sign(f(trial))
	if (length(unique(trial))==1) {
		warning('Seems Like f Might Not Have A Root In The Specified Interval.')
	}
}

if (precise) {
	n <- mpfr(1:(times+2), bits)
} else {
	n <- 1:(times+2)
}
n[1] <- guess1
n[2] <- guess2

for (i in 3:(times+2)){                                                       # The main event
	if (tol>0) {
			if ( abs(f(n[i-1])-f(n[i-2])) < tol | n[i-1]==666) {
				n[i] <- 666
			} else {
				n[i] <- n[i-1] - f(n[i-1])*(n[i-1]-n[i-2]) / (f(n[i-1])-f(n[i-2]))	
		} 
	} else {
		n[i] <- n[i-1] - f(n[i-1])*(n[i-1]-n[i-2]) / (f(n[i-1])-f(n[i-2]))	
	}
}
n <- n[1:(times+2)]
if (tol>0) {
	if ( any(abs(f(n))<tol) ) {
		for(i in min(which(abs(f(n))<tol)):times) {                             # Correct in case of accidental success
			n[i] <- n[min(which(abs(f(n))<tol))]; 
		}
	}
}
if (any(is.na(n))) {
	'We Had NaN(s) Somewhere. See'
	if (min(which(!is.na(n)))>4) {
		print(results[(min(which(is.na(n)))-4):min(which(is.na(n))),])
	} else {
		print(results[1:min(which(is.na(n))),])
	}
}
#
#	what to do with the results
#
if (precise) {
	results <- n
} else {
	results <- data.frame(Iteration=c(0,0:times), Estimated_Root=n)
}
if (!is.null(store)) {                                                        # Store results in a data frame inthe gloabl environment
	assign(store, results, envir = .GlobalEnv)
}
if (is.null(print) & !any(n==666)) {                                          # Show results in a table
	if (times<16) {
		print(results)
	} else {   
		print(tail(results))
	}
}
if (add) {
if (plot & !any(n==666)) {                                                    # Show results in a graph
	points(c(0,0:times), n, las=1, pch=c(7,7,rep(4,times=times)),
			xlab='Iterations', 
			ylab='Estimated Root', 
			main='Results of Secant Method')
	}
} else {
if (plot & !any(n==666)) {
	x11()                                                                   # Show results in a graph
	plot(c(0,0:times), n, las=1, pch=c(0,0,rep(1,times=times)),
			xlab='Iterations', 
			ylab='Estimated Root', 
			main='Results of Secant Method')
	}
}
if (any(n==666)) {                                                            # a half-assed attempt to disarm on a foot gun
	if (min(which(n==666))>4) {
		print(results[(min(which(n==666))-4):min(which(n==666)),])
	} else {
		print(results[1:min(which(n==666)),])
	}
	print(paste('We had f(a_n)=f(b_n) at some point ( tolerance', tol, '), indicated by 666'))
	}
if (is.logical(print)){
	if (print) results
}
#
#	end
#
}




#
##
###	that one dude's method
##
#

Newtons.Method <- function(f, f1, guess, times=100, store=NULL, plot=TRUE, print=NULL, check=TRUE, add=FALSE, tol=NULL, correct=F, precise=F, bits=250) {
#
#	warnings
#
if (!is.character(store)) {
	if (!is.null(store)) {stop("store Argument Must Be a Character String or NULL")}
}
if (!is.function(f)) {stop('f Must Be A Function Of x')}
if (!is.function(f1)) {stop('df Must Be A Function Of x')}
if (!is.null(tol)) {
	if (!is.numeric(tol)) {stop('tol Must Be Positive Real (and should be very small)')}
	if (tol<=0) {stop('tol Must Be Positive Real (and should be very small)')}
} else {
	if (precise) {tol <- 1e-10} else {tol <- 1e-30}
}
if (!is.numeric(guess)) {stop("Please Make First Initial Guess a Real Number")}
if (!is.null(dim(guess))) {stop("Please Make First Initial Guess a Real Number (not, e.g., a matrix)")}
if (length(guess)>1) {warning("guess1 Has Length Greater Than 1, Only The First Element Will Be Used")}

if (!is.numeric(times)) {stop("times Must Be A Positive Number")}
if (!is.null(dim(times))) {stop("times Must Be A Positive Number (not, e.g., a matrix)")}
if (times<1) {stop("times Must Be A Positive Number")}
if (length(times)>1) {warning("times Has Length Greater Than 1, Only The First Element Will Be Used")}
if (!floor(times[1])==times[1]) {warning("times Is Non-Integer, Will Be Rounded Down")}

if (!is.numeric(bits)) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (!is.null(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (is.list(dim(bits))) stop('bits Should Be Simply A Positive Integer (not, e.g., a matrix).')
if (length(bits)>1) warning('Only The First Entry Of bits Will Be Used.')
if (bits<1) stop('bits Must Be At Least 1 (and should probably be at least 53).')
if (!is.logical(precise)) stop('precise Must Be Logical.')	
if (precise & !'Rmpfr'%in%(.packages())) stop('In Order To Utilize The Setting precise=T, The Package Rmpfr Is Required.')
#
#	the actual algorithm
#
if (check) {
	ey0 <- guess - .3
	bee0 <- guess + .3
	trial <- seq(ey0, bee0, by=0.001)
	trial <- sign(f(trial))
	if (length(unique(trial))==1) {
		warning('Seems Like f Might Not Have A Root In The Specified Interval.')
	}
}
if (precise) {
	n <- mpfr(guess, bits)
} else {
	n <- guess
}
for(i in 2:(times+1)) {   
	if (is.na(f1(n[i-1]))) {
		n[i] <- NaN
	} else {
		if (abs(f1(n[i-1]))<tol | n[i-1]==666) { 
			n[i] <- 666
		} else {      
 			n[i] <- n[i-1] - f(n[i-1])/f1(n[i-1])
		}
	}
}
n <- n[1:(times+1)]
if (correct) {
	if ( any(abs(f(n[which(!is.na(n))]))<tol) ) {                                  # Correct in case of accidental success
		n <- n[which(!is.na(n))]
		for(i in min( which( abs(f(n))<tol )):(times+1)) {
			n[i] <- n[min(which(abs(f(n))<tol))]; 
    		}
	print( paste('Success. f(', n[min( which( abs(f(n[which(!is.na(n))]))<tol ) )], ')=0 to within tolerance', tol) )
	}
}                     
#
#	what to do with the results
#
if (precise) {
	results <- n
} else {
	results <- data.frame(Iteration=0:times, Estimated_Root=n)
}
if (!is.null(store)) {                                                        # Store results in a data frame inthe gloabl environment
	assign(store, results, envir = .GlobalEnv)
}
if (!(any(n==666) | any(is.na(n)))) {
if (is.null(print)) {                                          # Show results in a table
	if (times<16) {
		print(results)
	} else {   
		print(tail(results))
	}
}
if (add) {
if (plot) {                                                    # Show results in a graph
	points(0:times, n, las=1, pch=c(7,rep(4,times=times)),
			xlab='Iterations', 
			ylab='Estimated Root', 
			main='Results of Newtons Method')
	}
} else {
if (plot) {                                                    # Show results in a graph
	x11()
	plot(0:times, n, las=1, pch=c(0,rep(1,times=times)),
			xlab='Iterations', 
			ylab='Estimated Root', 
			main='Results of Newtons Method')
	}
}
}
if (any(is.na(n))) {
	print('We Had NaN(s) Somewhere. See')
	if (min(which(is.na(n)))>4) {
		print(results[(min(which(is.na(n)))-4):min(which(is.na(n))),])
	} else {
		print(results[1:min(which(is.na(n))),])
	}
} else {
	if (any(n==666)) {                                                            # a half-assed attempt to disarm on a foot gun
	if (min(which(n==666))>4) {
		print(results[(min(which(n==666))-4):min(which(n==666)),])
	} else {
		print(results[1:min(which(n==666)),])
	}
	print(paste('We had f1(p_n)=0 ( tolerance', tol,') at some point, indicated by 666'))
	}
}
if (is.logical(print)){
	if (print) results
}
#
#	end
#
}

f <- function(x) cos(x^cos(x)) - sin(x)
f1 <- function(x) -sin(x^cos(x))*x^cos(x)*(cos(x)/x - sin(x)*log(x)) - cos(x)
Newtons.Method(f, f1, 1, 100)




#
##
###	aitken's delta squared method
##
#

Aitkens <- function(n) {
	if (!is.numeric(n)) {stop('n Should Be A Numeric Vector (not, e.g., a matrix).')}
	if (!is.null(dim(n))) {stop('n Should Be A Numeric Vector (not, e.g., a matrix).')}
	if (length(n)<3) {stop('n Must Be A Numeric Vector Of Length 3 Or Longer.')}

	nhat <- 0; for (i in 1:(length(n)-2)) nhat[i] <- n[i]+n[i+2]-2*n[i+1]
	if ( any(nhat[!is.na(nhat)]==0) ) {
		print(paste('Cannot divide by zero. See place i =', min(which(nhat==0)), 'of n. Recall: we must divide by n_i + n_(i+2)-2n_(i+1)' ))
	} else {
		for (i in 1:(length(n)-2)) nhat[i] <- n[i]-(n[i]-n[i+1])^2/(n[i]+n[i+2]-2*n[i+1]) }
	first <- nhat[1]
	data.frame(Pre_Data=n[1:length(n)], Post_Data=c(first, first, nhat))

}





cff <- call.from.frame <- function(x, sto, label=TRUE){
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