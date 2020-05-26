%% Takes input information on binary trials (stimulus levels and correct data)
% outputs the t statistic and corresponding p value
% estimated values for mu and sigma also need to be inputted for fitting

function [t_ind,p_ind] = compareThresholds(X_sham,cor_sham,X_best,cor_best, mu_sham, mu_best, sigma_guess, N)


musi_bi_sham = jackknife(@just794_bi,X_sham,cor_sham,mu_sham,sigma_guess);
x = fminsearch(@(x) two_int_fit_simp(x, X_sham, cor_sham), [mu_sham, sigma_guess]);
mu_bi_sham = x(1);
sigma_bi_sham = x(2);
mu_se_bi_sham = sqrt( (N-1)/N*sum( (mu_bi_sham - musi_bi_sham(:,1)).^2) );
mu_jackData = musi_bi_sham(:,1);
mu_indDOF_sham=sum((mu_bi_sham - musi_bi_sham(:,1)).^2)^2/sum((mu_bi_sham - musi_bi_sham(:,1)).^4);
sigma_se_bi_sham = sqrt( (N-1)/N*sum( (sigma_bi_sham - musi_bi_sham(:,2)).^2 ) );

musi_bi_best = jackknife(@just794_bi,X_best,cor_best,mu_best,sigma_guess);
x = fminsearch(@(x) two_int_fit_simp(x, X_best, cor_best), [mu_best, sigma_guess]);
mu_bi_best = x(1);
sigma_bi_best = x(2);
mu_se_bi_best = sqrt( (N-1)/N*sum( (mu_bi_best - musi_bi_best(:,1)).^2) );
mu_jackData = musi_bi_best(:,1);
mu_indDOF_best=sum((mu_bi_best - musi_bi_best(:,1)).^2)^2/sum((mu_bi_best - musi_bi_best(:,1)).^4);
sigma_se_bi_best = sqrt( (N-1)/N*sum( (sigma_bi_best - musi_bi_best(:,2)).^2 ) );


mu_compDOF=(mu_se_bi_sham^2+mu_se_bi_best^2)^2/(1/mu_indDOF_sham*mu_se_bi_sham^4+1/mu_indDOF_best*mu_se_bi_best^4); % real DOF
t_ind=(mu_bi_sham-mu_bi_best)/sqrt(mu_se_bi_sham^2+mu_se_bi_best^2);
p_ind=tcdf(t_ind,mu_compDOF);
