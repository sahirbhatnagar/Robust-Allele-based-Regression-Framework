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
RA_indep <- function(g, y, z=NULL, HWE=F){
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
		if(class(z) == 'numeric'){
			if(var(z) <= 0){
				stop('Covariates are identical')
			}
		}else{
			temp_z_var <- apply(z, 2, var)
			if(any(temp_z_var <= 0)){
				stop('At least one of the covariates are identical')
			}
		}
	}
	if(var(g) <= 0){
		stop('All genotypes are identical')
	}
	if(var(y) <=0){
		stop('All phenotypes are identical')
	}
	gRA <- as.numeric(RA_split(g)); n <- length(y)
	yRA <- rep(y, each=2)
	p <- mean(gRA)
	num <- sum(yRA*(gRA-p))
	if(HWE){
		print('assuming HWE')
		sigma_g <- p*(1-p)
	}else{
		p2 <- mean(g==2)
		sigma_g <- p2-p^2+p*(1-p)
	}
	sigma_y <- sum(yRA^2)-2*n*mean(yRA)^2
	return(num^2/sigma_g/sigma_y)
}

RA_indep <- function(g, y, z){
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
	return(RA_stat[1,1])
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

### RA test for tri-allelic locus



### RA HWE test for siblings 



