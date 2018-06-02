#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


NumericVector cbin(NumericVector x, int n, int nbin, double min, double max){
  NumericVector output(nbin);
  double skip = (nbin-1.0)/(max-min);
  for(int i=0; i<n; i++) output[floor((x[i]-min)*skip)] +=1;
  return output;
}


NumericVector cbin2(NumericVector x, int n, int nbin, double min, double max){
  NumericVector output(nbin);
  double skip = (nbin-1.0)/(max-min);
  int fl = 0;
  for(int i=0; i<n; i++){
    fl = floor((x[i]-min)*skip);
    if(fl<(nbin-1)){
      output[fl+1] += (x[i]-min)*skip-fl;
      output[fl] += fl+1+(min-x[i])*skip;
    }
    else output[nbin-1] +=1;
  } 
  return output;
}

// [[Rcpp::export]]

NumericVector f_kde_mix(NumericVector x, double h, int n, int ord, NumericVector betas){
  NumericVector srt = clone(x);
  std::sort(srt.begin(), srt.end());
  NumericVector output(n);
  double denom;
  NumericMatrix L(ord + 1, n);
  NumericMatrix R(ord + 1, n);
  for(int i=0; i<=ord; i++) L(i,0) = pow(-srt[0], i);
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = pow(-srt[i],j) + exp((srt[i-1]-srt[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((srt[n-i-1]-srt[n-i])/h)*(pow(srt[n-i],j)+R(j,n-i));
    }
  }
  for(int orddo=0; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = n*pow(h, orddo+1);
    for(int i=0; i<n; i++){
      for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*(pow(srt[i], orddo-j)*L(j,i)+pow(-srt[i],orddo-j)*R(j,i))/denom;
    }
  }
  return output;
}


// [[Rcpp::export]]

NumericVector f_kde_mix_grid(NumericVector x, double h, int n, int ord, NumericVector betas, int ngrid, double mn, double MX){
  NumericVector srt = clone(x);
  std::sort(srt.begin(), srt.end());
  NumericMatrix L(ord + 1, n);
  NumericMatrix R(ord + 1, n);
  for(int i=0; i<=ord; i++) L(i,0) = pow(-srt[0], i);
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = pow(-srt[i],j) + exp((srt[i-1]-srt[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((srt[n-i-1]-srt[n-i])/h)*(pow(srt[n-i],j)+R(j,n-i));
    }
  }
  NumericVector nx = cbin(srt, n, ngrid, mn, MX);
  NumericVector output(ngrid);
  double skip = (MX-mn)/(ngrid-1.0);
  double xe;
  double exp_mult;
  double denom;
  int count;
  for(int orddo=0; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = n*pow(h, orddo+1);
    count = round(nx[0]);
    for(int i=0; i<ngrid; i++){
      xe = mn+i*skip;
      if(count==0){
        exp_mult = exp((xe-srt[0])/h);
        output[i] += betas[orddo]*pow(srt[0]-xe, orddo)/denom*exp_mult;
        for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*pow(-xe,orddo-j)*R(j,0)/denom*exp_mult;
      }
      else{
        exp_mult = exp((srt[count-1]-xe)/h);
        for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*(pow(xe, orddo-j)*L(j,count-1)*exp_mult+pow(-xe,orddo-j)*R(j,count-1)/exp_mult)/denom;
      }
      count += round(nx[i]);
    }
  }
  return output;
}

// [[Rcpp::export]]


NumericVector f_kde_mix_binned(NumericVector x, double h, int n, int ord, NumericVector betas, int nbin, double min, double max){
  double skip = (max-min)/(nbin-1.0);
  NumericVector nx = cbin2(x, n, nbin, min+skip/1000000, max+skip/1000000);
  double eskip = exp(-skip/h);
  NumericMatrix L(ord + 1, nbin);
  NumericMatrix R(ord + 1, nbin);
  double denom=0;
  for(int i=1; i<nbin; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = eskip*L(j,i-1) + nx[i]*pow(-(min+i*skip),j);
      R(j,nbin-i-1) = eskip*(nx[nbin-i]*pow(min+(nbin-i)*skip,j) + R(j,nbin-i));
    }
  }
  NumericVector output(nbin);
  double add = 0;
  double tot = 0;
  for(int orddo=0; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = pow(h, orddo+1);
    for(int i=0; i<nbin; i++){
      add = 0;
      for(int j=0; j<=orddo; j++){
        add += betas[orddo]/denom*coefs[j]*(pow(min+i*skip, orddo-j)*L(j,i)+pow(-min-i*skip,orddo-j)*R(j,i));
      }
      output[i] += add;
      tot += add;
    }
  }
  for(int i=0; i<nbin; i++) output[i] /= tot*skip;
  return output;
}


// [[Rcpp::export]]

NumericVector df_kde(NumericVector x, double h, int n, int ord){
  double alfact = 1;
  for(int i=2; i<=(ord+1); i++){
    alfact *= i;
  }
  NumericVector srt = clone(x);
  std::sort(srt.begin(), srt.end());
  NumericVector output(n);
  double denom = 2.0*alfact*pow(h,ord+2)*n;
  NumericMatrix L(ord + 1, n);
  NumericMatrix R(ord + 1, n);
  for(int i=0; i<=ord; i++) L(i,0) = pow(-srt[0], i);
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = pow(-srt[i],j) + exp((srt[i-1]-srt[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((srt[n-i-1]-srt[n-i])/h)*(pow(srt[n-i],j)+R(j,n-i));
    }
  }
  NumericVector coefs(ord + 1);
  coefs[0] = coefs[ord] = 1;
  if(ord>1){
    double num = 1;
    for(int j=2; j<=ord; j++) num *= j;
    double denom1 = 1;
    double denom2 = num/ord;
    for(int i=2; i<=ord; i++){
      coefs[i-1] = num/denom1/denom2;
      denom1 *= i;
      denom2 /= (ord-i+1);
    }
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<=ord; j++) output[i] -= coefs[j]*(pow(srt[i], ord-j)*L(j,i)-pow(-srt[i],ord-j)*R(j,i))/denom;
  }
  return output;
}
