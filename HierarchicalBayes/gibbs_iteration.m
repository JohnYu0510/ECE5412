function output_param = gibbs_iteration(input_param,data,sim_param)
%--------------------------------------------------------------------------
% One step of Gibbs iterations to update the parameter estimators
% input_param: struct, parameters simulated from the last iteration
% data: struct, market data of stock returns and factor returns
% sim_param: struct, simulation parameters
%--------------------------------------------------------------------------
% Generate samples of factor moments
output_param.mu_f = mean(data.f_t)+sqrt(input_param.Omega_f/sim_param.T)*randn;
output_param.Omega_f = 1/wishrnd(2/(sum((data.f_t-input_param.mu_f).^2)),sim_param.T/2);
% Generate samples of factor loadings (alphas)
y_jt_star = data.y_jt - data.f_t*input_param.beta';
term1 = (input_param.theta_0*data.zj_alpha./input_param.Lambda_0+sim_param.T*mean(y_jt_star)'./input_param.v);
term2 = 1./(1./input_param.Lambda_0+sim_param.T./input_param.v);
output_param.alpha = term1.*term2+sqrt(term2).*randn(sim_param.J,1);
% Generate samples of factor loadings (betas)
y_jt_star = data.y_jt - repmat(input_param.alpha',sim_param.T,1);
term1 = input_param.theta_1*data.zj_beta./input_param.Lambda_1+(y_jt_star'*data.f_t)./input_param.v;
term2 = 1./(1/input_param.Lambda_1+sum(data.f_t.^2)./input_param.v);
output_param.beta = term1.*term2+sqrt(term2).*randn(sim_param.J,1);
% Generate samples of factor loadings (taus)
imbedded_mh = false;
Sj_over_T = mean((y_jt_star-data.f_t*input_param.beta').^2)';
term1 = sim_param.T/2*log(Sj_over_T)+input_param.psi*data.zj_tau/input_param.delta;
term2 = 1/(sim_param.T/2+1/input_param.delta);
output_param.tau = nan(sim_param.J,1);
mh_start = term1*term2+sqrt(term2)*randn(sim_param.J,1);
mh_mean = term1*term2; mh_var = term2;
if imbedded_mh
    for j = 1:sim_param.J
        mh_pdf = @(x)exp(-sim_param.T/2*x-(sim_param.T/2*(Sj_over_T(j))).*exp(-x)...
            -0.5/input_param.delta*(x-input_param.psi*data.zj_tau(j))^2);
        mh_proppdf = @(x,y)normpdf(x,mh_mean(j),sqrt(mh_var));
        mh_proprnd = @(x)mh_mean(j)+sqrt(mh_var)*randn;
        mh_samples = mhsample(mh_start(j),100,'pdf',mh_pdf,'proppdf',mh_proppdf,'proprnd',mh_proprnd);
        output_param.tau(j) = mh_samples(end);
    end
else
    output_param.tau = mh_start;
end
output_param.v = exp(output_param.tau);
% Generate samples of hyperparameters
output_param.theta_0 = regress(input_param.alpha,data.zj_alpha)+sqrt(input_param.Lambda_0/(data.zj_alpha'*data.zj_alpha))*randn;
output_param.theta_1 = regress(input_param.beta,data.zj_beta)+sqrt(input_param.Lambda_1/(data.zj_beta'*data.zj_beta))*randn;
output_param.psi = regress(input_param.tau,data.zj_tau)+sqrt(input_param.delta/(data.zj_tau'*data.zj_tau))*randn;
a0 = 1; a1 = 0.1;
output_param.Lambda_0 = 1/gamrnd(a0+sim_param.J/2,a1+sum((input_param.alpha-input_param.theta_0*data.zj_alpha).^2)/2);
output_param.Lambda_1 = 1/gamrnd(a0+sim_param.J/2,a1+sum((input_param.beta-input_param.theta_1*data.zj_beta).^2)/2);
output_param.delta = 1/gamrnd(a0+sim_param.J/2,a1+sum((input_param.tau-input_param.psi*data.zj_tau).^2)/2);
% Compute asset moments
output_param.E_y = input_param.alpha + input_param.beta*input_param.mu_f;
output_param.Var_y = input_param.beta'*input_param.Omega_f*input_param.beta+diag(input_param.v);
end