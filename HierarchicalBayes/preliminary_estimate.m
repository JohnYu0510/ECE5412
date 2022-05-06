function initial = preliminary_estimate(data,sim_param)
%--------------------------------------------------------------------------
% Obtain the preliminary estimates from least squares estimation
% data: struct, market data of stock returns and factor returns
% sim_param: struct, simulation parameters
%--------------------------------------------------------------------------
initial.alpha = nan(sim_param.J,1); initial.beta = nan(sim_param.J,1); initial.v = nan(sim_param.J,1);
for j = 1:sim_param.J
    coeff = regress(data.y_jt(:,j),data.regressor);
    initial.alpha(j) = coeff(1); initial.beta(j) = coeff(2);
    initial.v(j) = var(data.y_jt(:,j)-data.regressor*coeff);
end
initial.tau = log(initial.v);
[initial.theta_0,~,res] = regress(initial.alpha,data.zj_alpha);
initial.Lambda_0 = var(res);
[initial.theta_1,~,res] = regress(initial.beta,data.zj_beta);
initial.Lambda_1 = var(res);
[initial.psi,~,res] = regress(initial.tau,data.zj_tau);
initial.delta = var(res);
initial.mu_f = mean(data.f_t); initial.Omega_f = var(data.f_t);
initial.E_y = initial.alpha+initial.mu_f*initial.beta;
initial.Var_y = initial.beta'*initial.Omega_f*initial.beta+diag(initial.v);
end