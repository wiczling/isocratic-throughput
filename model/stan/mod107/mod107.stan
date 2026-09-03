functions {
// credit: http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/

vector lower_tri(matrix mat) {
  int d = rows(mat);
  int n = d * (d - 1) / 2;
  vector[n] out;
  int idx = 1;

  for (r in 2:d) {
    for (c in 1:(r - 1)) {
      out[idx] = mat[r, c];
      idx += 1;
    }
  }

  return out;
}
  
  real lkj_corr_cholesky_point_lower_tri2_lpdf(matrix cor_L, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(cor[2,1]  | point_mu_lower, point_scale_lower);
    return(lpdf);
 }
 
  real lkj_corr_cholesky_point_lower_tri_lpdf(matrix cor_L, vector point_mu_lower, vector point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(lower_tri(cor) | point_mu_lower, point_scale_lower);
    return(lpdf);
 }
 
  real funlogki_fi(real logkw, real S1, real S2, real gamma, real alpha, real pHmpKa, real xlogkw, real xS1, real xS2, real fi) {
  // isocratic Neue model
  // hinge function for dewetting
   real log10 = 2.302585092994046;
   real logfix1  = -(1+xS2)*xS1*fi/(1+xS2*fi);
   real logfix2  = -(1+S2)*S1 *fi/(1+S2*fi);
   real t2 = (pHmpKa  + alpha * fi) * log10;

   real a = (xlogkw + gamma*fi + logfix1) * log10;
   real c = (logkw + logfix2)*log10 + t2;
   
   real logk = log_sum_exp(a, c)/log10-log1p_exp(t2)/log10;
   //return logkd - delta*log(1+exp((logkd-logk)/delta));
   return logk;
  }
  
real funlogki(real logkw, real S1, real logS2, real gamma, real fid, real alpha, real pHmpKa, real dlogkw, real dS1, real dlogS2, real fi){
    real S2  = pow(10, logS2);
    real xS1 = S1+dS1;
    real xlogkw = logkw+dlogkw;
    real xS2  = pow(10, logS2+dlogS2);  
    
    real logkd     = funlogki_fi(logkw, S1, S2, gamma, alpha, pHmpKa, xlogkw, xS1, xS2, fid);
    real logk_raw  = funlogki_fi(logkw, S1, S2, gamma, alpha, pHmpKa, xlogkw, xS1, xS2,  fi);
    real w = inv_logit((fid - fi) / 0.001);
    return (w*logkd+(1-w)*logk_raw);
   }
  
  real partial_sum(array[] int ind, int start, int end, 
                   vector logkobs,
                   array[] int analyte,
                   array[] int column,
                   array[] real fi,
                   array[,] vector param,
                  real sigma) {
                     
    real lp = 0;

    for (z in start : end) {
      
    int i = analyte[z];
    int c = column[z];
    real y_hat = funlogki(param[i,c,1], param[i,c,2], param[i,c,3], param[i,c,4], param[i,c,5], 
                          param[i,c,6], param[i,c,7], param[i,c,8], param[i,c,9], param[i,c,10], fi[z]);
    lp = lp + normal_lpdf(logkobs[z] | y_hat,  sigma);
    }
    return lp;
  }
}

data {
  int nAnalytes; 
  int nColumns;
  int nObs;                 
  array[nObs] int analyte;
  array[nObs] int column;
  array[nObs] real fi; 
  array[nAnalytes] real logP;
  array[nColumns] real logCL;
  array[nAnalytes] real dissociated;
  array[nAnalytes] real pHmpKa;
  vector[nObs] logkobs;                 // observed log retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  int grainsize = 1;
  array[nObs] int ind = rep_array(1, nObs);

}

parameters {
  
  //typical values
  real logkw_hat;                   // typical logkw
  real S1_hat;                      // effect of ACN on logkw
  real logS2_hat;                   // typical value of logS2
  vector[2] beta_a_hat;             // effect of logP on logkw, S1, logS2
  vector[3] beta_c_hat;             // effect of carbon load on logkw, S1, logS2
  real gamma_hat;
  real alpha_hat;
  real fid_hat;
  
  real dlogS2_hat;
  real dlogkw_hat;
  real dS1_hat;
  
  // analyte  
  vector<lower=0>[3] omega_a; 
  cholesky_factor_corr[3] L_a;
  //column
  vector<lower=0>[8] omega_c;
  cholesky_factor_corr[3] L_c;
  //analyte x column  

  vector[3] m_logomega_d;
  vector<lower=0>[3] sd_logomega_d;
  cholesky_factor_corr[3] L_d;
  
  real<lower=0> sigma;
  matrix[3, nColumns] logomega_dStd;
  matrix[3, nAnalytes] etaStd;
  matrix[8, nColumns] kappaStd;
  
  matrix[nAnalytes,nColumns] dzeta1;
  matrix[nAnalytes,nColumns] dzeta2;
  matrix[nAnalytes,nColumns] dzeta3;
}

transformed parameters {
   array[nAnalytes,nColumns] vector[10] param;
   array[nAnalytes] vector[10] itheta;
   array[nColumns] vector[3] ctheta;
   matrix[nColumns, 3] logomega_d;
   matrix[nAnalytes, 3] eta;
   matrix[nColumns, 8] kappa;
   vector[nColumns] omega_d1; 
   vector[nColumns] omega_d2;
   vector[nColumns] omega_d3;
   
   logomega_d =  diag_pre_multiply(sd_logomega_d, L_d * logomega_dStd)'; 
   
   omega_d1 = exp(m_logomega_d[1]+logomega_d[,1]);
   omega_d2 = exp(m_logomega_d[2]+logomega_d[,2]);
   omega_d3 = exp(m_logomega_d[3]+logomega_d[,3]);
  
   eta =  diag_pre_multiply(omega_a, L_a * etaStd)'; 
   
   kappa[,1:3]=  diag_pre_multiply(omega_c[1:3], L_c * kappaStd[1:3,])'; 
   kappa[,4]=  omega_c[4] * kappaStd[4,]';
   kappa[,5]=  omega_c[5] * kappaStd[5,]';
   kappa[,6]=  omega_c[6] * kappaStd[6,]';
   kappa[,7]=  omega_c[7] * kappaStd[7,]';
   kappa[,8]=  omega_c[8] * kappaStd[8,]';

   
  for (i in 1:nAnalytes) {
   itheta[i,1] = logkw_hat + beta_a_hat[1]*(logP[i]-2.5);
   itheta[i,2] = S1_hat    + beta_a_hat[2]*(logP[i]-2.5);
   itheta[i,3] = logS2_hat ;
   itheta[i,4] = gamma_hat;
   itheta[i,5] = fid_hat;
   itheta[i,6] = alpha_hat;
   itheta[i,7] = pHmpKa[i];
   itheta[i,8] = dlogkw_hat;
   itheta[i,9] = dS1_hat;
   itheta[i,10] = dlogS2_hat;
   }
    
   for (c in 1:nColumns) {
   ctheta[c,1] = beta_c_hat[1]*(logCL[c]-1);
   ctheta[c,2] = beta_c_hat[2]*(logCL[c]-1);
   ctheta[c,3] = beta_c_hat[3]*(logCL[c]-1);
   }
   
  for (i in 1:nAnalytes) {
   for (c in 1:nColumns) {
    param[i,c,1] = itheta[i,1] + ctheta[c,1] + eta[i,1] + kappa[c,1] +  omega_d1[c]*dzeta1[i,c];
    param[i,c,2] = itheta[i,2] + ctheta[c,2] + eta[i,2] + kappa[c,2] +  omega_d2[c]*dzeta2[i,c]; 
    param[i,c,3] = itheta[i,3] + ctheta[c,3] + eta[i,3] + kappa[c,3] +  omega_d3[c]*dzeta3[i,c];
    param[i,c,4] = itheta[i,4] + kappa[c,4]; //gamma
    param[i,c,5] = itheta[i,5] + kappa[c,5]; //fid
    param[i,c,6] = itheta[i,6]; //alpha (pH vs fi)
    param[i,c,7] = itheta[i,7]; //pH-pKa
    param[i,c,8]  = itheta[i,8] + kappa[c,6];   //dlogkw_hat
    param[i,c,9]  = itheta[i,9] + kappa[c,7];   //dS1_hat
    param[i,c,10] = itheta[i,10] + kappa[c,8];  //dlogS2_hat
   }}
}

model {
   logkw_hat ~ normal(2.5, 1);
   S1_hat    ~ normal(5, 2);
   logS2_hat ~ normal(0.30103, 0.05);
   fid_hat   ~ normal(0,0.025);
   beta_a_hat[{1}] ~ normal(0.75, 0.25);
   beta_a_hat[{2}] ~ normal(0.5, 0.25);
   beta_c_hat ~ normal(0.5, 0.25);
   gamma_hat  ~ normal(1, 0.5);
   alpha_hat  ~ normal(3, 0.25);
   dlogkw_hat ~ normal(-1, 0.5);
   dS1_hat ~ normal(0, 0.5);
   dlogS2_hat ~ normal(0, 0.5);
   omega_a[{1,2}] ~ normal(0, 1);   //logkw and S1
   omega_a[{3}]   ~ normal(0, 0.1); //logS2
   omega_c[{1,2}] ~ normal(0, 1);   //logkw and S1
   omega_c[{3}]   ~ normal(0, 0.1); //logS2 
   omega_c[{4}]   ~ normal(0, 0.5);
   omega_c[{5}]   ~ normal(0, 0.025);
   
   omega_c[{6}]   ~ normal(0, 0.1);
   omega_c[{7}]   ~ normal(0, 0.1);
   omega_c[{8}]   ~ normal(0, 0.1);
   
   to_vector(logomega_dStd) ~ normal(0, 1);
   
   m_logomega_d ~ normal(-0.693,1); 
   sd_logomega_d ~ normal(0,1);

  L_a ~ lkj_corr_cholesky_point_lower_tri([0.75, 0.00, 0.00]', [0.25,0.25,0.25]');
  L_c ~ lkj_corr_cholesky_point_lower_tri([0.75, 0.75, 0.75]', [0.25,0.25,0.25]');
  L_d ~ lkj_corr_cholesky(5);
  
  to_vector(etaStd) ~ normal(0, 1);
  to_vector(kappaStd) ~ normal(0, 1);

  to_vector(dzeta1) ~ normal(0,1); 
  to_vector(dzeta2) ~ normal(0,1);
  to_vector(dzeta3) ~ normal(0,1);
   
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
   // logkobs ~ student_t(7,logkx, sigma);
     target += reduce_sum(partial_sum, ind, grainsize, logkobs,
                          analyte, column, fi, param,
                          sigma);
  }
}

generated quantities {
  vector[nObs] logkx;
  corr_matrix[3] rho_a;
  corr_matrix[3] rho_c;
  corr_matrix[3] rho_d;
  
  rho_a = multiply_lower_tri_self_transpose(L_a);
  rho_c = multiply_lower_tri_self_transpose(L_c);
  rho_d = multiply_lower_tri_self_transpose(L_d);
  
   for (z in 1 : nObs) {
    int i = analyte[z];
    int c = column[z];
    logkx[z] = funlogki(param[i,c,1], param[i,c,2], param[i,c,3], param[i,c,4], param[i,c,5], param[i,c,6],
                        param[i,c,7], param[i,c,8], param[i,c,9], param[i,c,10], fi[z]);
  }
}
