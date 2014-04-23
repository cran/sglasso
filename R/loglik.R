loglik <- function(object, size = 2){
	if(missing(size)) stop("sample size is not specified")
	if(size <= 0) stop("size must be a positive integer")
	if(size != ceiling(size)) stop("size must be an integer value")
	p <- dim(object$S)[1]
	rho <- object$rho
	nv <- object$nv
	ne <- object$ne
	w <- object$w
	th_e <- object$theta[(nv + 1):(nv + ne), , drop = FALSE]
	l1_th_e <- apply(th_e, 2, function(x) {
					 id <- which(abs(x) > 0)
					 sum(w[id] * abs(x[id]))
					 }
					 )
	out <- Kh(object)
	out <- unlist(lapply(out, det))
	k <- 0.5 * size
	out <- k * (log(out) - p + rho * l1_th_e)
	out
}
