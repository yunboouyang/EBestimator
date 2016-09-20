library("EbayesThresh")
library("wavethresh")
library("Rsolnp")

set.seed(100)
#Data generation process
num=200;

x0=array(0,dim=c(6,num,1000))
b=c(2,3,4,5,6,7);
for (time in 1:num){
  n=1; #diffenrent situations
  for (i in 1:6){
    truemu=c(rep(0,900), rep(b[i],90),rep(10,10));
    x0[n,time,]=rnorm(1000, truemu);
    n=n+1;
  }
}



#EBMW
VB<-function(x,alpha,sigma,kappa,print=TRUE){ #x is a vector,  alpha, sigma and kappa are parameters
  n=length(x);
  E0=rep(0,n);
  cut=1.96*sqrt(1/kappa);  #initial cutoff values
  E=ifelse(abs(x)>cut,0,1);
  num=0;
  while(max(abs(E-E0))>10^-3){
    num=num+1;
    E0=E;
    alphaw=alpha*n+sum(E);
    betaw=n-sum(E)+1;
    E=exp(log(sqrt(1+kappa*sigma^2))+digamma(alphaw)-digamma(betaw)-kappa/2*x^2)/(exp(log(sqrt(1+kappa*sigma^2))+digamma(alphaw)-digamma(betaw)-kappa/2*x^2)+1)
  }
  estimate=(1-E)*x;
  check=pnorm(abs(x)/sqrt(sigma^2/(1+kappa*sigma^2)))/(1+exp(log(sqrt(1+kappa*sigma^2))+digamma(alphaw)-digamma(betaw)-kappa/2*x^2));
  select=ifelse(check>0.5,1,0)
  result=data.frame(prob=E,estimate=estimate,select=select);
  return(result);
}

dis.CD = function(x, prior.mass, h=1){
  # prior on mu is a discrete mixture over values provided in prior.mass
  # bandwith = h
  n=length(prior.mass)
  tmp = outer(x, prior.mass, '-'); 
  tmp = exp(-tmp^2/(2*h)); 
  tmp = tmp/rowSums(tmp); 
  return(tmp %*% prior.mass)
  
}





sparseVBDP=function(x,alpha, sigma, w, T0=10){
  n=length(x);
  lambda=1/sigma^2;
  phi0=matrix(rep(0,n*T0),nrow=n,ncol=T0);
  phi=matrix(rep(1/T0,n*T0),nrow=n,ncol=T0);
  while(max(abs(phi-phi0))>10^-3){
    phi0=phi;
    gamma1=1+colSums(phi)[1:(T0-1)];
    gamma2=alpha+rev(cumsum(rev(colSums(phi)))[1:(T0-1)]);
    tau1=(t(phi)%*%x)[1:T0,];
    tau2=lambda+colSums(phi);
    odds=log(w)-log(1-w)+log(sqrt(1/lambda*colSums(phi)+1))-tau1^2/(2*tau2);
    p=exp(odds)/(1+exp(odds))
    d1=c(digamma(gamma1)-digamma(gamma1+gamma2),0);
    d2=digamma(gamma2)-digamma(gamma1+gamma2);
    d3=d1+c(0,cumsum(d2));
    d=d3-0.5*(1-p)*((tau1/tau2)^2+1/tau2);
    S=outer(x,(1-p)*tau1/tau2,'*')+outer(rep(1,n),d,'*');
    E=exp(S);
    phi=E/rowSums(E);
  }
  mean=c(0,tau1/tau2);
  
  newphi=cbind(phi%*%p,phi%*%diag(1-p));
  number=max.col(newphi);
  prior=mean[number];
  csize=length(unique(number));
  return(list(prior=prior,csize=csize));
}




opt=function(f){   #This is the objective function for EBKW
  g=A%*%f;
  return(-sum(log(g)))
}


dis.CD2 = function(x, prior, grid, h=1){  #This function is used to construct EBKW
  # every grid has a probability mass
  # bandwith = h
  tmp = outer(x, grid, '-'); 
  tmp = exp(-tmp^2/(2*h))%*%diag(prior); 
  tmp = tmp/rowSums(tmp); 
  return(as.vector(tmp %*% grid))
  
}








MSE=array(0,dim=c(10,6,num));
MAE=array(0,dim=c(10,6,num));

series=seq(1,1000,1);
u=seq(from=-11.92,to=12,by=0.08)
thresh=sqrt(2*log(1000))


for (time in 1:num){
  n=1; #diffenrent situations
  for (i in 1:6){
      x=x0[n,time,];
      #EBMW
      truemu=c(rep(0,900), rep(b[i],90),rep(10,10));
      mu=VB(x, 0.10, 10,0.99);  
      MSE[1,n,time]=sum((mu$estimate-truemu)^2);
      MAE[1,n,time]=sum(abs(mu$estimate-truemu));
      
      
      #EBKW
      A=exp(-0.5*outer(x,u,'-')^2)/(2*pi);
      prior=solnp(pars=rep(1/300,300),fun=opt,eqfun=sum,eqB=1,LB=rep(0,300),control=list(trace=0))$pars
      mu=dis.CD2(x,prior,u)
      MSE[2,n,time]=sum((mu-truemu)^2);
      MAE[2,n,time]=sum(abs(mu-truemu));
      
      #Johnstone and Silverman
      mu=ebayesthresh(x, a=NA,prior="laplace", sdev=1);
      MSE[3,n,time]=sum((mu-truemu)^2);
      MAE[3,n,time]=sum(abs(mu-truemu));
      
      #SURE
      suth=sure(x);
      mu=threshld(x, suth, hard=F); 
      MSE[4,n,time]=sum((mu-truemu)^2);
      MAE[4,n,time]=sum(abs(mu-truemu));
      
      #Soft thresholding;
      mu=threshld(x, thresh, hard=F); 
      MSE[5,n,time]=sum((mu-truemu)^2);
      MAE[5,n,time]=sum(abs(mu-truemu));
      
      #Hard thresholding;
      mu=threshld(x, thresh, hard=T);   
      MSE[6,n,time]=sum((mu-truemu)^2);
      MAE[6,n,time]=sum(abs(mu-truemu));
      
      
      #FDR
      q=0.01;
      seq=qnorm(1-q/2*series/1000);
      y=sort(abs(x),decreasing=T);
      binary=ifelse(y>=seq,1,0);
      m=max(binary*series);
      if (m==0) {mu=rep(0,1000);}
      else {mu=threshld(x, seq[m], hard=T);}   #Hard thresholding;
      MSE[7,n,time]=sum((mu-truemu)^2);
      MAE[7,n,time]=sum(abs(mu-truemu));
      
      #FDR
      q=0.1;
      seq=qnorm(1-q/2*series/1000);
      y=sort(abs(x),decreasing=T);
      binary=ifelse(y>=seq,1,0);
      m=max(binary*series);
      if (m==0) {mu=rep(0,1000);}
      else {mu=threshld(x, seq[m], hard=T);}   #Hard thresholding;
      MSE[8,n,time]=sum((mu-truemu)^2);
      MAE[8,n,time]=sum(abs(mu-truemu));
      
      
      #FDR
      q=0.4;
      seq=qnorm(1-q/2*series/1000);
      y=sort(abs(x),decreasing=T);
      binary=ifelse(y>=seq,1,0);
      m=max(binary*series);
      if (m==0) {mu=rep(0,1000);}
      else {mu=threshld(x, seq[m], hard=T);}   #Hard thresholding;
      MSE[9,n,time]=sum((mu-truemu)^2);
      MAE[9,n,time]=sum(abs(mu-truemu));
      
      
      #VBDP bayes factor
      results=sparseVBDP(x,1,6,0.01)
      prior=results$prior;
      mu=dis.CD(x,prior);
      
      MSE[10,n,time]=sum((mu-truemu)^2);
      MAE[10,n,time]=sum(abs(mu-truemu));
      
      
      
      n=n+1;
  }
  
  print(time)
  
}

MSE=round(apply(MSE,c(1,2),mean))
MAE=round(apply(MAE,c(1,2),mean))

save(MSE,MAE,file="/home/youyang4/estimation code/1000derror.RData")


