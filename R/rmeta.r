rmeta <- function(y,v,model="RE",gamma=0.01){

	alpha0 <- seq(gamma,(1-gamma),by=gamma)

	R <- NULL

	if(model=="RE"){

	for(i in 1:length(alpha0)){

		alpha <- alpha0[i]
		
		rma1 <- rma(y, v, method="DL")

		mu1 <- as.numeric(rma1$beta)
		V1 <- as.numeric(rma1$tau2.f)

		L1 <- function(mu){

			Vi <- (2*pi*(v+V1))^(-0.5)
			Ai <- -alpha*(y-mu)^2
			Bi <- 2*(v+V1)
			si <- Vi^alpha * exp(Ai/Bi)

			-sum( si/alpha - Vi^alpha/((1+alpha)^1.5) )	# minus pseudo-likelihood

		}

		L2 <- function(V){

			Vi <- (2*pi*(v+V))^(-0.5)
			Ai <- -alpha*(y-mu2)^2
			Bi <- 2*(v+V)
			si <- Vi^alpha * exp(Ai/Bi)

			-sum( si/alpha - Vi^alpha/((1+alpha)^1.5) )	# minus pseudo-likelihood

		}

		repeat{

			mu2 <- optimize(L1, c(-5,5))$minimum
			V2 <- optimize(L2, c(0.01,5))$minimum

			rmse <- (abs(mu1-mu2) + abs(V1-V2))/2
			if(rmse<10^-4) break
		
			mu1 <- mu2; V1 <- V2

		}
	
		###
	
		Ci <- v*v/(v+V2)
		Vi <- (2*pi*(v+V2))^(-0.5)
		Di <- Vi^(2*alpha) / ((2*alpha+1)^1.5)
		Ei <- 2*Vi^alpha / ((alpha+1)^1.5)

		g1 <- (v*V2)/(v+V2)
		g2 <- Ci * (Di - Ei + 1)
	
		# Exc <- sum(g2)/sum(g1)		# criterion for the selection of alpha by Sugasawa (2020)
		Exc <- sum(g1+g2)		# another criterion for the selection of alpha; simply summing up the mean of squared predicted errors 

		###

		Fi <- 2*(alpha*(y-mu2)^2 - (v+V2)) / ((v+V2)*(v+V2))
		Gi <- dnorm(y, mean=mu2, sd=sqrt(v+V2))^alpha
		Hi <- (y-mu2)^2 / ((v+V2)*(v+V2))
		Ii <- dnorm(y, mean=mu2, sd=sqrt(v+V2))^(2*alpha)
		Hn <- sum(Fi*Gi + Hi*Ii)

		###

		Jb <- sum( Vi^alpha/(v+V2) ) / ( length(y) * ((alpha+1)^1.5) )
		Kb <- sum( Vi^(2*alpha)/(v+V2) ) / ( length(y) * ((2*alpha+1)^1.5) )
		Vb <- Kb/(length(y)*Jb*Jb)
		SEb <- sqrt(Vb)

		###

		R <- rbind(R, c(alpha,mu2,V2,SEb,Hn,Exc))

	}
	
	w1 <- which.min(R[,5])
	R1 <- R[w1,]
	
	CI <- R1[2] + c(-1,1)*qnorm(.975)*R1[4]
	Z <- -abs(R1[2]/R1[4])
	P <- 2*pnorm(Z)
	
	###

	Vi <- (2*pi*(v+R1[3]))^(-0.5)
	Ai <- -R1[1]*(y-R1[2])^2
	Bi <- 2*(v+R1[3])
	hi <- Vi^R1[1] * exp(Ai/Bi)
			
	Ci <- hi/((v+R1[3])^2)
	Di <- R1[3] + v - R1[1]*(y-R1[2])*(y-R1[2])
	Ei <- Ci*Di
	
	wi <- Ei/sum(Ei)		# approximate contribution rate
	
	V3 <- rma(y, v, method="REML")$tau2
	ui <- 1/(v+V3); ui <- ui/sum(ui)	# contribution rate of REML estimate
	
	# xi <- wi/ui		# relative change of the contribution rates

	study <- 1:length(wi)

	W <- data.frame(study,ui,wi) #,xi)

	###
	
	R2 <- list(mu=R1[2],se=R1[4],CI=CI,P=P,tau2=R1[3],alpha=R1[1],W=W)
	
	return(R2)
	
	}
	
	if(model=="FE"){

	for(i in 1:length(alpha0)){

		alpha <- alpha0[i]

		rma1 <- rma(y, v, method="EE")

		mu1 <- as.numeric(rma1$beta)

		L1 <- function(mu){

			Vi <- (2*pi*v)^(-0.5)
			Ai <- -alpha*(y-mu)^2
			Bi <- 2*v
			si <- Vi^alpha * exp(Ai/Bi)

			-sum( si/alpha - Vi^alpha/((1+alpha)^1.5) )	# minus pseudo-likelihood

		}

		mu2 <- optimize(L1, c(-5,5))$minimum
	
		###
	
		Fi <- 2*(alpha*(y-mu2)^2 - v) / (v*v)
		Gi <- dnorm(y, mean=mu2, sd=sqrt(v))^alpha
		Hi <- (y-mu2)^2 / (v*v)
		Ii <- dnorm(y, mean=mu2, sd=sqrt(v))^(2*alpha)
		Hn <- sum(Fi*Gi + Hi*Ii)

		###

		Vi <- (2*pi*v)^(-0.5)
		Jb <- sum( Vi^alpha/v ) / ( length(y) * ((alpha+1)^1.5) )
		Kb <- sum( Vi^(2*alpha)/v ) / ( length(y) * ((2*alpha+1)^1.5) )
		Vb <- Kb/(length(y)*Jb*Jb)
		SEb <- sqrt(Vb)

		###

		R <- rbind(R, c(alpha,mu2,SEb,Hn))

	}
	
	w1 <- which.min(R[,4])
	R1 <- R[w1,]
	
	CI <- R1[2] + c(-1,1)*qnorm(.975)*R1[3]
	Z <- -abs(R1[2]/R1[3])
	P <- 2*pnorm(Z)
	
	###
	
	###

	Vi <- (2*pi*v)^(-0.5)
	Ai <- -R1[1]*(y-R1[2])^2
	Bi <- 2*v
	hi <- Vi^R1[1] * exp(Ai/Bi)
			
	Ci <- hi/(v^2)
	Di <- v - R1[1]*(y-R1[2])*(y-R1[2])
	Ei <- Ci*Di
	
	wi <- Ei/sum(Ei)		# approximate contribution rate
	
	ui <- 1/v; ui <- ui/sum(ui)	# contribution rate of REML estimate
	
	# xi <- wi/ui		# relative change of the contribution rates

	study <- 1:length(wi)

	W <- data.frame(study,ui,wi) #,xi)
	
	###
	
	R2 <- list(mu=R1[2],se=R1[3],CI=CI,P=P,alpha=R1[1],W=W)
	
	return(R2)
	
	}
	
}

