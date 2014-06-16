summary.sglasso <- function (object, size, k = c("bic","aic"), digits = max(3, getOption("digits") - 3), ...){
	if(missing(size)) stop("size is not specified")
	if(size <= 0) stop("k must be a positive integer")
	if(is.character(k)){
		type <- match.arg(k)
		k <- ifelse(type == "bic", log(size), 2)
	} else {
		if(k < 0) stop("k must be greater than zero")
		type <- "GoF"
	}
	rho <- object$rho
	df <- object$df
	action <- make_action(object$theta)
	ll <- loglik(object, size)
	gof <- - 2 * ll + k * df
	gof_rank <- rank(gof)
	gof_mark <- vector(mode = "character", length = object$nrho)
	gof_min <- which(gof_rank == 1)
	gof_mark[gof_min] <- "<-"
	theta_gof <- object$theta[, gof_min]
	theta_gof <- theta_gof[abs(theta_gof) > 0]
	tbl <- data.frame(action, rho, ll, df, gof, gof_rank, gof_mark)
	names(tbl) <- c("Sequence", "rho", "log-lik", "df", type, "rank", "")
	tbl.format <- format(tbl, digits = digits)
	cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	cat("\n==============================================\n")
	cat("\nCoefficients estimated by", type, "criterion ( k =", k, "):\n\n")
	print.default(format(theta_gof, digits = digits), print.gap = 2, quote = FALSE, ...)
	cat("\n", type, ": ", min(tbl.format[, 5]))
	cat("\n\nAlgorithm", object$algorithm, "with exit =", object$conv, "( n. iterations =", object$n, ")", "\n\n")
	out <- list(table = tbl, theta_gof = theta_gof)
	invisible(out)
}
