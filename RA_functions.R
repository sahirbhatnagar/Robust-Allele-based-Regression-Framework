### partition genotypes to alleles
RA_split <- function(g){
	sapply(g, function(x) ((x==1)*sample(c(0,1), size=2) + (x != 1)*c(x/2, x/2)))
}

### Test of HWE using RA framework
RA_HWE <- function(g){
	if(any(is.na(g))){
		print(paste('Removing', sum(is.na(g)), 'NA observations'))
		g <- g[!is.na(g)]
	}
	if(length(g) <=1){
		stop('Fewer than two non-NA observations')
	}
	if(!all(g==1 | g==0 | g==2)){
		stop('Genotypes contain values other than 0, 1 or 2')
	}
	if(var(g) <= 0){
		stop('All genotypes are identical')
	}
	n <- length(g); n2 <- sum(g==2)
	p2 <- n2/n; p <- mean(g)/2
	rho <- (p2-p^2)/p/(1-p)
	pchisq(n*rho^2, df=1, lower.tail=F)
}


### Association test for independent samples
RA_assoc_indep <- function(g, y, z=NULL, HWE=F){
	## 1: screen for NAs, check var(G), var(Y), var(Z)
	if(any(is.na(g))){
		stop('Genotypes contain NAs')
	}
	if(any(is.na(y))){
		stop('Phenotypes contain NAs')
	}
	if(!is.null(z)){
		if(any(is.na(z))){
			stop('Covariates contain NAs')
		}
	}else{
		if(class(z)=='numeric'){
			z <- matrix(z, ncol=1)
		}
		z_var <- apply(z, 2, var)
		if(any(z_var <= 0)){
			stop('variance of Z is 0')
		}
	}
	if(var(g) <= 0){
		stop('variance of G is 0')
	}
	if(class(y) == 'numeric'){
		y <- matrix(y, ncol=1)
	}
	y_var <- apply(y, 2, var)
	if(any(y_var <=0)){
		stop('variance of Y is 0')
	}
	if(!is.null(z)){
		return(RA_indep_Z(g=g, y=y, z=z))
	}else{
		return(RA_indep(g=g, y=y))	
	}
}

RA_indep_Z <- function(g, y, z){
	n <- length(g)
	g_RA <- RA_split(g); gbar <- mean(g_RA)
	g_vec <- matrix(g_RA-gbar, ncol=1)
	one_vec <- matrix(rep(1, length(g_RA)), ncol=1)
	y_RA <- apply(y, 2, function(x) rep(x, each=2))
	z_RA <- apply(z, 2, function(x) rep(x, each=2))
	X <- cbind(one_vec, y_RA, z_RA)
	z_bar <- apply(z_RA, 2, mean)
	z_bar <- matrix(z_bar, ncol=1)
	var_gamma <- t(z_RA)%*%z_RA - z_bar%*%t(z_bar)
	gamma_hat <- solve(var_gamma)%*%t(z_RA)%*%g_vec
	score_fn <- t(y_RA)%*%(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)
	score_vec <- matrix(c(0, score_fn[,1], rep(0, ncol(z))), ncol=1)
	sigma_hat <- t(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)%*%(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)/length(g_RA)
	rho_mat <- (g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)[,1]
	rho_mat <- matrix(rho_mat, byrow=T, ncol=2)
	rho_hat <- sum(rho_mat[,1]*rho_mat[,2])/n/sigma_hat
	RA_stat <- t(score_vec)%*%solve(t(X)%*%X)%*%score_vec/sigma_hat/(1+rho_hat)
	RA_pval <- pchisq(RA_stat[1,1], df=nrow(score_fn), lower.tail=F)
	return(RA_pval)
}

RA_assoc_indep <- function(g, y, z=NULL, HWE=F){
	## 1: screen for NAs, check var(G), var(Y), var(Z)
	if(any(is.na(g))){
		stop('Genotypes contain NAs')
	}
	if(any(is.na(y))){
		stop('Phenotypes contain NAs')
	}
	if(!is.null(z)){
		if(any(is.na(z))){
			stop('Covariates contain NAs')
		}
	}else{
		if(class(z)=='numeric'){
			z <- matrix(z, ncol=1)
		}
		z_var <- apply(z, 2, var)
		if(any(z_var <= 0)){
			stop('variance of Z is 0')
		}
	}
	if(var(g) <= 0){
		stop('variance of G is 0')
	}
	if(class(y) == 'numeric'){
		y <- matrix(y, ncol=1)
	}
	y_var <- apply(y, 2, var)
	if(any(y_var <=0)){
		stop('variance of Y is 0')
	}
	if(!is.null(z)){
		return(RA_indep_Z(g=g, y=y, z=z))
	}else{
		return(RA_indep(g=g, y=y))	
	}
}

RA_indep_Z <- function(g, y, z){
	n <- length(g)
	g_RA <- RA_split(g); gbar <- mean(g_RA)
	g_vec <- matrix(g_RA-gbar, ncol=1)
	one_vec <- matrix(rep(1, length(g_RA)), ncol=1)
	y_RA <- apply(y, 2, function(x) rep(x, each=2))
	z_RA <- apply(z, 2, function(x) rep(x, each=2))
	X <- cbind(one_vec, y_RA, z_RA)
	z_bar <- apply(z_RA, 2, mean)
	z_bar <- matrix(z_bar, ncol=1)
	var_gamma <- t(z_RA)%*%z_RA - z_bar%*%t(z_bar)
	gamma_hat <- solve(var_gamma)%*%t(z_RA)%*%g_vec
	score_fn <- t(y_RA)%*%(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)
	score_vec <- matrix(c(0, score_fn[,1], rep(0, ncol(z))), ncol=1)
	sigma_hat <- t(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)%*%(g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)/length(g_RA)
	rho_mat <- (g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)[,1]
	rho_mat <- matrix(rho_mat, byrow=T, ncol=2)
	rho_hat <- sum(rho_mat[,1]*rho_mat[,2])/n/sigma_hat
	RA_stat <- t(score_vec)%*%solve(t(X)%*%X)%*%score_vec/sigma_hat/(1+rho_hat)
	RA_pval <- pchisq(RA_stat[1,1], df=nrow(score_fn), lower.tail=F)
	return(RA_pval)
}

RA_indep <- function(g, y){
	n <- length(g)
	g_RA <- RA_split(g); gbar <- mean(g_RA)
	g_vec <- matrix(g_RA-gbar, ncol=1)
	y_RA <- apply(y, 2, function(x) rep(x, each=2))
	score_fn <- apply(y_RA, 2, function(x) sum(x*g_vec))
	score_vec <- matrix(c(0, score_fn[,1])), ncol=1)
	sigma_hat <- gbar*(1-gbar)
	g_mat <- matrix(g_vec, ncol=2, byrow=T)
	rho_hat <- sum(g_mat[,1]*g_mat[,2])/n/sigma_hat
	one_vec <- matrix(rep(1, n*2), ncol=1)
	y_mat <- as.matrix(cbind(one_vec, y_RA))
	info_inv <- solve(t(y_mat)%*%y_mat)
	RA_stat <- t(score_vec)%*%info_inv%*%score_vec/sigma_hat/(1+rho_hat)
	RA_pval <- pchisq(RA_stat[1,1], df=nrow(score_fn), lower.tail=F)
	return(RA_pval)
}

### RA test of siblings + Independent

RA_pedi_assoc <- function(gIndep, gSib, yMat_indep, yMat_sib, phi=0.25, thres=10^(-8)){
	para_vec <- fast_rho(gIndep=gIndep, gSib=gSib, phi=phi, crit=thres)
	alpha <- para_vec[1]; sigma <- para_vec[2]; rho <- para_vec[3]
	score_fn <- sapply(1:ncol(yMat_indep), function(k) mix_score_fn(gIndep=gIndep, gSib=gSib, yIndep=yMat_indep[,k], ySib=yMat_sib[,k], alpha=alpha, sigma=sigma, rho=rho))
	fisher_inverse <- multi_info(yMat_indep=yMat_indep, yMat_sib=yMat_sib, sigma=sigma, rho=rho)
	score_fn <- matrix(c(0, score_fn), ncol=1)
	test_stat <- t(score_fn)%*%fisher_inverse%*%score_fn
	return(test_stat[1,1])
}

