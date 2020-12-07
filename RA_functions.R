library(Matrix)

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
		return(RA_indep_Z(g=g, y=y, z=z, HWE=HWE))
	}else{
		return(RA_indep(g=g, y=y, HWE=HWE))	
	}
}

RA_indep_Z <- function(g, y, z, HWE=F){
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
	if(HWE){
		rho_hat <- 0
	}else{
		rho_mat <- (g_vec - (z_RA - one_vec%*%t(z_bar))%*%gamma_hat)[,1]
	rho_mat <- matrix(rho_mat, byrow=T, ncol=2)
	rho_hat <- sum(rho_mat[,1]*rho_mat[,2])/n/sigma_hat
	}
	RA_stat <- t(score_vec)%*%solve(t(X)%*%X)%*%score_vec/sigma_hat/(1+rho_hat)
	RA_pval <- pchisq(RA_stat[1,1], df=nrow(score_fn), lower.tail=F)
	return(RA_pval)
}

RA_indep <- function(g, y, HWE=F){
	n <- length(g)
	g_RA <- RA_split(g); gbar <- mean(g_RA)
	g_vec <- matrix(g_RA-gbar, ncol=1)
	y_RA <- apply(y, 2, function(x) rep(x, each=2))
	score_fn <- apply(y_RA, 2, function(x) sum(x*g_vec))
	score_vec <- matrix(c(0, score_fn), ncol=1)
	sigma_hat <- gbar*(1-gbar)
	g_mat <- matrix(g_vec, ncol=2, byrow=T)
	if(HWE){
		rho_hat <- 0
	}else{
		rho_hat <- sum(g_mat[,1]*g_mat[,2])/n/sigma_hat
	}
	one_vec <- matrix(rep(1, n*2), ncol=1)
	y_mat <- as.matrix(cbind(one_vec, y_RA))
	info_inv <- solve(t(y_mat)%*%y_mat)
	RA_stat <- t(score_vec)%*%info_inv%*%score_vec/sigma_hat/(1+rho_hat)
	RA_pval <- pchisq(RA_stat[1,1], df=length(score_fn), lower.tail=F)
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
	pval <- pchisq(test_stat[1,1], df=length(score_fn), lower.tail=F)
	return(pval)
}

fast_alpha <- function(gIndep, gSib, rho, phi=0.25){
	a <- phi*(1+rho); f <- length(gSib); n <- length(gIndep)
	alpha <- (sum(gIndep)/(1+rho) + sum(gSib)/(1+2*a))/(2*n/(1+rho) + 2*f/(1+2*a))
	return(alpha)
}

fast_sigma <- function(gIndep, gSib, rho, phi=0.25){
	alpha <- fast_alpha(gIndep=gIndep, gSib=gSib, rho=rho, phi=phi)
	a <- phi*(1+rho); f <- length(gSib); n <- length(gIndep)
	RA_ind <- RA_split(gIndep)-alpha
	RA_sib <- RA_split(gSib)-alpha
	ind_mat <- matrix(RA_ind, ncol=2, byrow=T)
	ind_kernel <- sum(RA_ind^2)-2*rho*sum(ind_mat[,1]*ind_mat[,2])
	sib_mat <- matrix(RA_sib, ncol=4, byrow=T)
	sib_kernel <- (1-2*a^2)*sum(RA_sib^2) + 4*a^2*sum(sib_mat[,1]*sib_mat[,2] + sib_mat[,3]*sib_mat[,4]) - 2*a*sum((sib_mat[,1]+sib_mat[,2])*(sib_mat[,3]+sib_mat[,4]))
	sigma <- (ind_kernel/(1-rho^2) + sib_kernel/(1-4*a^2))/(2*n + 2*f)
	return(c(alpha, sigma))
}

rho_second <- function(gIndep, gSib, rho, phi=0.25){
	a <- phi*(1+rho); f <- length(gSib); n <- length(gIndep)
	para <- fast_sigma(gIndep, gSib, rho=rho)
	alpha <- para[1]; sigma <- para[2]
	RA_ind <- RA_split(gIndep)-alpha
	ind_mat <- matrix(RA_ind, ncol=2, byrow=T)
	sib_mat <- matrix((gSib-2*alpha), ncol=2, byrow=T)
	ind_kernel <- n*rho/(1-rho^2) - (rho*sum(RA_ind^2)-(1+rho^2)*sum(ind_mat[,1]*ind_mat[,2]))/sigma/(1-rho^2)^2
	sib_kernel <- 2*f*a/(1-4*a^2) - (2*a*sum(gSib^2) - (1+4*a^2)*sum(sib_mat[,1]*sib_mat[,2]))/sigma/(1-4*a^2)^2
	rho_score <- ind_kernel + sib_kernel/4
	ind_2nd <- n*(1+rho^2)/(1-rho^2)^2 - ((1+3*rho^2)*sum(RA_ind^2)-2*rho*(3+rho^2)*sum(ind_mat[,1]*ind_mat[,2]))/sigma/(1-rho^2)^3
	sib_2nd <- 2*f*(1+4*a^2)/(1-4*a^2)^2 - (2*(1+12*a^2)*sum(gSib^2) - 8*a*(3+4*a^2)*sum(sib_mat[,1]*sib_mat[,2]))/sigma/(1-4*a^2)^3
	rho_2nd <- ind_2nd + sib_2nd/16
	return(c(rho_score, rho_2nd, alpha, sigma))
}

fast_rho <- function(gIndep, gSib, phi=0.25, crit=10^(-8)){
	if(sum(gIndep, gSib) <= 0){
		return(NA)
	}
	rho_temp <- 0
	temp_values <- rho_second(gIndep=gIndep, gSib=gSib, rho=rho_temp, phi=phi)
	i <- 0
	while(abs(temp_values[1]) > crit){
		rho_temp <- rho_temp - temp_values[1]/temp_values[2]
		temp_values <- rho_second(gIndep=gIndep, gSib=gSib, rho=rho_temp, phi=phi)
	}
	return(c(temp_values[3:4], rho_temp))
}

mix_score_fn <- function(gIndep, gSib, yIndep, ySib, alpha, rho, sigma, phi=0.25){
	a <- phi*(1+rho)
	score_indep <- sum(yIndep*(gIndep-2*alpha))/sigma/(1+rho)
	y_sib_mat <- matrix(ySib, ncol=2, byrow=T)
	g_sib_mat <- matrix(gSib-2*alpha, ncol=2, byrow=T)
	score_sib <- (sum(ySib*(gSib-2*alpha)) - 2*a*sum(y_sib_mat[,1]*g_sib_mat[,2] + y_sib_mat[,2]*g_sib_mat[,1]))/sigma/(1-4*a^2)
	score_total <- score_indep + score_sib
	return(score_total)
}


multi_info <- function(yMat_indep, yMat_sib, sigma, rho, phi=0.25){
	k <- ncol(yMat_indep)
	n <- nrow(yMat_indep); f <- nrow(yMat_sib); a <- phi*(1+rho)
	info_mat <- matrix(ncol=(k+1), nrow=(k+1))
	info_mat[1,1] <- 2*n/(1+rho) + 2*f/(1+2*a)
	info_mat[1,2:(k+1)] <- sapply(1:k, function(x) 2*sum(yMat_indep[,x])/(1+rho) + 2*sum(yMat_sib[,x])/(1+2*a))
	for(i in 1:k){
		for(l in i:k){
			info_mat[(i+1),(l+1)] <- mix_info_vec(u_indep=yMat_indep[,i], u_sib=yMat_sib[,i], w_indep=yMat_indep[,l], w_sib=yMat_sib[,l], rho=rho, a=a)
		}
	}
	final_info <- forceSymmetric(info_mat)
	final_info <- final_info/sigma
	return(solve(final_info))
}

mix_info_vec <- function(u_indep, u_sib, w_indep, w_sib, rho, a ){
	u_sib_mat <- matrix(u_sib, ncol=2, byrow=T)
	w_sib_mat <- matrix(w_sib, ncol=2, byrow=T)
	indep_sum <- 2*sum(u_indep*w_indep)/(1+rho)
	sib_sum <- 2*(sum(u_sib*w_sib) -2*a*(sum(u_sib_mat[,1]*w_sib_mat[,2]) + sum(u_sib_mat[,2]*w_sib_mat[,1])))/(1-4*a^2)
	return(indep_sum + sib_sum)
}
