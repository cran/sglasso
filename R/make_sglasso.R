make_sglasso <- function(object, call, algorithm, S, mask){
	if(object$conv == 1){
		nr <- object$nrho
		theta <- object$th[, 1:nr, drop = FALSE]
		grd <- object$grd[, 1:nr, drop = FALSE]
		df <- object$df[1:nr]
		rho <- object$rho[1:nr]
	}	else {
		theta <- object$th
		grd <- object$grd
		df <- object$df
		rho <- object$rho
	}
	w <- object$w	
	nv <- object$nv
	ne <- object$ne
	nrho <- ifelse(is.null(object$nrho), 1, object$nrho)
	obj <- list(call = call, nv = nv, ne = ne, theta = theta, w = w, df = df, rho = rho, grd = grd, 
				nstep = object$nstep, nrho = nrho, algorithm = algorithm, tol = object$tol, S = S, 
				mask = mask, n = object$n, conv = object$conv)
	class(obj) <- "sglasso"
	obj
}
