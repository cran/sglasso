sglasso <- function(S, mask, w = NULL, min_rho = 1.0e-02, nrho = 50, nstep = 1.0e+05, algorithm = c("ccd","ccm"), tol = 1.0e-05){
	this.call <- match.call()
	if(class(S) == "matrix") S <- as(S, "dspMatrix")
	if(is.null(dimnames(S)[[1]])) dimnames(S)[[1]] <- paste("X",1:dim(S)[1],sep="")
	if(is.null(dimnames(S)[[2]])) dimnames(S)[[2]] <- dimnames(S)[[1]]
	if(missing(mask)) stop("mask is not specified. See the documentation for more details")
	if(any(dim(S) != dim(mask))) stop("dim(S) is different from dim(mask)")
	if(storage.mode(mask) != "character") storage.mode(mask) <- "character"
	if(any(is.na(mask[upper.tri(mask)]))) mask[is.na(mask)] <- "zero"
	if(min_rho < 0) stop("min_rho can not be a negative value. See the documentation for more details")
	if(nrho < 0) stop("nrho can not be a negative value. See the documentation for more details")
	if(nstep < 0) stop("nstep can not be a negative value. See the documentation for more details")
	algo <- match.arg(algorithm)
	if(tol < 0) stop("tol can not be a negative value. See the documentation for more details")
	out.fit <- sglasso.fit(Sv = S@x, mask = mask, w = w, nrho = nrho, min_rho = min_rho, nstep = nstep, algorithm = algo, tol = tol)
	out.fit <- make_sglasso(object = out.fit, call = this.call, algo, S = S, mask = mask)
	out.fit
}
