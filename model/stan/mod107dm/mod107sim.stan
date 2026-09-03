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
  array[nAnalytes] real pHmpKa;
  vector[nObs] logkobs;                 // observed log retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
  
  int nObs_sim ;
  int nAnalytes_sim;
  vector[nObs_sim] logkobs_sim; 
  array[nObs_sim] real fi_sim;
  array[nAnalytes_sim] real logP_sim;
  array[nAnalytes_sim] real pHmpKa_sim;
  array[nObs_sim] int analyte_sim;
  array[nObs_sim] int column_sim;
}

transformed data{
  int grainsize = 1;
  array[nObs_sim] int ind = rep_array(1, nObs_sim);
  
    
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
  
  vector<lower=0>[3] omega_a; 
  matrix[3,3] L_a;
  vector<lower=0>[8] omega_c;
  matrix[3,3] L_c;
  vector[3] m_logomega_d;
  vector<lower=0>[3] sd_logomega_d;
  matrix[3,3] L_d;
  real<lower=0> sigma;
  
  matrix[8, nColumns] kappaStd;
  matrix[3, nColumns] logomega_dStd;
  
logkw_hat=3.16117193435;
S1_hat=3.948875903925;
logS2_hat=0.32338299042;
fid_hat=0.0120711463189716;
alpha_hat=4.732614796125;
dlogkw_hat=-2.3795255041;
dS1_hat=2.159345496275;
dlogS2_hat=-0.8619173294475;
gamma_hat=3.064047275775;
beta_a_hat[1]=0.8307343083275;
beta_a_hat[2]=0.7400127906975;
beta_c_hat[1]=0.736271500472602;
beta_c_hat[2]=0.543006545386475;
beta_c_hat[3]=0.280537252878425;
omega_a[1]=0.9920858812275;
omega_a[2]=0.915700610125;
omega_a[3]=0.345018226315;
omega_c[1]=0.426174629735;
omega_c[2]=0.300703928135;
omega_c[3]=0.10105087003525;
omega_c[4]=1.264422655325;
omega_c[5]=0.011942497371825;
omega_c[6]=0.2106394137275;
omega_c[7]=0.465998469495;
omega_c[8]=0.16239544502175;
m_logomega_d[1]=-1.7850142501725;
m_logomega_d[2]=-2.006281283625;
m_logomega_d[3]=-2.7095769630925;
sd_logomega_d[1]=0.7126428989575;
sd_logomega_d[2]=0.63528297169;
sd_logomega_d[3]=1.43973800846;
L_a[1,1]=1;
L_a[2,1]=0.938189647125;
L_a[3,1]=0.444635396515;
L_a[1,2]=0;
L_a[2,2]=0.34374071537;
L_a[3,2]=-0.190673242265498;
L_a[1,3]=0;
L_a[2,3]=0;
L_a[3,3]=0.865593157845;
L_c[1,1]=1;
L_c[2,1]=0.82898312192;
L_c[3,1]=0.70891481467;
L_c[1,2]=0;
L_c[2,2]=0.53494352729;
L_c[3,2]=0.271308854510979;
L_c[1,3]=0;
L_c[2,3]=0;
L_c[3,3]=0.59057939116;
L_d[1,1]=1;
L_d[2,1]=0.430138685829688;
L_d[3,1]=0.594501766501475;
L_d[1,2]=0;
L_d[2,2]=0.8777093623125;
L_d[3,2]=0.206376319215216;
L_d[1,3]=0;
L_d[2,3]=0;
L_d[3,3]=0.7177368089375;
sigma=0.0207151967305;

kappaStd[1,1]=-1.28423180212848;
kappaStd[2,1]=1.54946771513085;
kappaStd[3,1]=-0.47781687671892;
kappaStd[4,1]=0.447258132270887;
kappaStd[5,1]=1.36077316054825;
kappaStd[6,1]=1.5815478679275;
kappaStd[7,1]=0.160228710726621;
kappaStd[8,1]=0.178277713799283;
kappaStd[1,2]=-2.6143841243325;
kappaStd[2,2]=-1.399397357136;
kappaStd[3,2]=-0.497094164800351;
kappaStd[4,2]=0.857046110795875;
kappaStd[5,2]=-0.379802302443603;
kappaStd[6,2]=1.13990570908101;
kappaStd[7,2]=-0.853295898906087;
kappaStd[8,2]=1.09382569273085;
kappaStd[1,3]=-0.112320407128718;
kappaStd[2,3]=-0.75066830527325;
kappaStd[3,3]=-0.299507505882089;
kappaStd[4,3]=0.43226459466219;
kappaStd[5,3]=-0.21062703607177;
kappaStd[6,3]=2.0538558348675;
kappaStd[7,3]=-0.956650315767583;
kappaStd[8,3]=1.6642834804035;
kappaStd[1,4]=0.384122446003902;
kappaStd[2,4]=0.531994871615821;
kappaStd[3,4]=1.9218087769951;
kappaStd[4,4]=1.05413953789025;
kappaStd[5,4]=-0.368164980694855;
kappaStd[6,4]=-0.663574082600905;
kappaStd[7,4]=-0.967981910303275;
kappaStd[8,4]=0.456439838259445;
kappaStd[1,5]=0.187206454349418;
kappaStd[2,5]=-0.944337578066625;
kappaStd[3,5]=0.190419211668597;
kappaStd[4,5]=-0.476303465111397;
kappaStd[5,5]=1.22509177173425;
kappaStd[6,5]=-0.746914204054165;
kappaStd[7,5]=1.29291099894761;
kappaStd[8,5]=-2.0433360845;
kappaStd[1,6]=0.182012121208097;
kappaStd[2,6]=0.102852486735128;
kappaStd[3,6]=-0.0495785455507425;
kappaStd[4,6]=0.311379858399413;
kappaStd[5,6]=-0.46730048176328;
kappaStd[6,6]=-0.267433199283563;
kappaStd[7,6]=-0.0727608784377643;
kappaStd[8,6]=-0.0952189468321125;
kappaStd[1,7]=-0.03039249578692;
kappaStd[2,7]=0.125292006264403;
kappaStd[3,7]=1.09977445383784;
kappaStd[4,7]=1.2763612440102;
kappaStd[5,7]=0.9418258480003;
kappaStd[6,7]=-0.74499326859846;
kappaStd[7,7]=2.3496055099825;
kappaStd[8,7]=-1.72213180725812;
kappaStd[1,8]=0.430694224702838;
kappaStd[2,8]=-0.283389245456202;
kappaStd[3,8]=-0.142785078154534;
kappaStd[4,8]=2.851907477075;
kappaStd[5,8]=-0.452644955589955;
kappaStd[6,8]=-1.46434033231425;
kappaStd[7,8]=2.7354970742675;
kappaStd[8,8]=-0.610366078313803;
kappaStd[1,9]=0.78102405066685;
kappaStd[2,9]=0.07279597872622;
kappaStd[3,9]=0.83226854808569;
kappaStd[4,9]=0.41188883636587;
kappaStd[5,9]=-0.381946315108783;
kappaStd[6,9]=0.295176130693922;
kappaStd[7,9]=-0.76460446567615;
kappaStd[8,9]=0.51512315206663;
kappaStd[1,10]=0.154303144334289;
kappaStd[2,10]=0.63742106564017;
kappaStd[3,10]=0.27742664841715;
kappaStd[4,10]=-0.347213811109708;
kappaStd[5,10]=-0.476868845966412;
kappaStd[6,10]=0.0160819647366325;
kappaStd[7,10]=-2.06219630206125;
kappaStd[8,10]=0.40886502698647;
kappaStd[1,11]=0.63933983683685;
kappaStd[2,11]=0.808925755173103;
kappaStd[3,11]=-1.09355084203003;
kappaStd[4,11]=1.130533334071;
kappaStd[5,11]=0.0251685932466882;
kappaStd[6,11]=-1.3182202666933;
kappaStd[7,11]=0.477135788436776;
kappaStd[8,11]=0.096056678362452;
kappaStd[1,12]=1.17562102424;
kappaStd[2,12]=-0.067525701687445;
kappaStd[3,12]=-0.468591596938015;
kappaStd[4,12]=0.981039061746475;
kappaStd[5,12]=-0.353970026129213;
kappaStd[6,12]=-0.262331017434763;
kappaStd[7,12]=2.033152980665;
kappaStd[8,12]=-0.411547341179885;
kappaStd[1,13]=1.12914862289675;
kappaStd[2,13]=-0.181788277310252;
kappaStd[3,13]=-0.794363369215025;
kappaStd[4,13]=1.12902504464178;
kappaStd[5,13]=-0.392987124011758;
kappaStd[6,13]=-0.774670821950428;
kappaStd[7,13]=0.634469658243084;
kappaStd[8,13]=-0.048790382590547;


logomega_dStd[1,1]=0.956767266591675;
logomega_dStd[2,1]=0.261291391598186;
logomega_dStd[3,1]=0.956484682986153;
logomega_dStd[1,2]=1.6435846827225;
logomega_dStd[2,2]=1.3462975642385;
logomega_dStd[3,2]=-0.238743783784091;
logomega_dStd[1,3]=0.39506475540913;
logomega_dStd[2,3]=-0.487425164117615;
logomega_dStd[3,3]=0.732234426158778;
logomega_dStd[1,4]=-0.38343822103908;
logomega_dStd[2,4]=0.558267819463602;
logomega_dStd[3,4]=-0.0702985943931025;
logomega_dStd[1,5]=0.757734104748575;
logomega_dStd[2,5]=1.04595132346142;
logomega_dStd[3,5]=-0.0726147876757025;
logomega_dStd[1,6]=-1.671451176715;
logomega_dStd[2,6]=-0.867917010550845;
logomega_dStd[3,6]=-1.05466656778922;
logomega_dStd[1,7]=0.0443727037609695;
logomega_dStd[2,7]=0.410273980955431;
logomega_dStd[3,7]=0.72140228074085;
logomega_dStd[1,8]=0.752730921541055;
logomega_dStd[2,8]=-0.522615092274022;
logomega_dStd[3,8]=-0.0367518832857675;
logomega_dStd[1,9]=-0.447578763290857;
logomega_dStd[2,9]=-1.68390497608725;
logomega_dStd[3,9]=-0.495160230608167;
logomega_dStd[1,10]=-1.14821366838657;
logomega_dStd[2,10]=-0.141033036850826;
logomega_dStd[3,10]=-1.16021125440705;
logomega_dStd[1,11]=-0.112706195462717;
logomega_dStd[2,11]=-0.698412751531792;
logomega_dStd[3,11]=-0.317713746701953;
logomega_dStd[1,12]=-1.6596546592;
logomega_dStd[2,12]=-0.211545838876337;
logomega_dStd[3,12]=-0.878737305114785;
logomega_dStd[1,13]=-2.0482473244925;
logomega_dStd[2,13]=-0.25542032168981;
logomega_dStd[3,13]=-0.048962051868835;

}

parameters {
  matrix[3, nAnalytes_sim] etaStd;
  matrix[nAnalytes_sim,nColumns] dzeta1;
  matrix[nAnalytes_sim,nColumns] dzeta2;
  matrix[nAnalytes_sim,nColumns] dzeta3;
}

transformed parameters {
   array[nAnalytes_sim,nColumns] vector[10] param;
   
   {
   array[nAnalytes_sim] vector[10] itheta;
   array[nColumns] vector[3] ctheta;
   matrix[nColumns, 3] logomega_d;
   matrix[nAnalytes_sim, 3] eta;
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

   
  for (i in 1:nAnalytes_sim) {
   itheta[i,1] = logkw_hat + beta_a_hat[1]*(logP_sim[i]-2.5);
   itheta[i,2] = S1_hat    + beta_a_hat[2]*(logP_sim[i]-2.5);
   itheta[i,3] = logS2_hat ;
   itheta[i,4] = gamma_hat;
   itheta[i,5] = fid_hat;
   itheta[i,6] = alpha_hat;
   itheta[i,7] = pHmpKa_sim[i];
   itheta[i,8] = dlogkw_hat;
   itheta[i,9] = dS1_hat;
   itheta[i,10] = dlogS2_hat;
   }
    
   for (c in 1:nColumns) {
   ctheta[c,1] = beta_c_hat[1]*(logCL[c]-1);
   ctheta[c,2] = beta_c_hat[2]*(logCL[c]-1);
   ctheta[c,3] = beta_c_hat[3]*(logCL[c]-1);
   }
   
  for (i in 1:nAnalytes_sim) {
   for (c in 1:nColumns) {
    param[i,c,1]  = itheta[i,1] + ctheta[c,1] + eta[i,1] + kappa[c,1] + omega_d1[c]*dzeta1[i,c];
    param[i,c,2]  = itheta[i,2] + ctheta[c,2] + eta[i,2] + kappa[c,2] + omega_d2[c]*dzeta2[i,c]; 
    param[i,c,3]  = itheta[i,3] + ctheta[c,3] + eta[i,3] + kappa[c,3] + omega_d3[c]*dzeta3[i,c];
    param[i,c,4]  = itheta[i,4] + kappa[c,4]; //gamma
    param[i,c,5]  = itheta[i,5] + kappa[c,5]; //fid
    param[i,c,6]  = itheta[i,6];              //alpha (pH vs fi)
    param[i,c,7]  = itheta[i,7];              //pH-pKa
    param[i,c,8]  = itheta[i,8] + kappa[c,6]; //dlogkw_hat
    param[i,c,9]  = itheta[i,9] + kappa[c,7]; //dS1_hat
    param[i,c,10] = itheta[i,10]+ kappa[c,8]; //dlogS2_hat
   }}
   
   }
}

model {
  to_vector(etaStd) ~ normal(0,1);
  to_vector(dzeta1) ~ normal(0,1); 
  to_vector(dzeta2) ~ normal(0,1);
  to_vector(dzeta3) ~ normal(0,1);
  
  if (run_estimation == 1) {
   // logkobs ~ student_t(7,logkx, sigma);
     target += reduce_sum(partial_sum, ind, grainsize, logkobs_sim,
                          analyte_sim, column_sim, fi_sim, param,
                          sigma);
  }
}
generated quantities {
}

