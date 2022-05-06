function [real_param,data] = simulate_data(data,sim_param)
%--------------------------------------------------------------------------
% Simulate factor loadings and returns data
% data: struct, observed z_j of companies
% sim_param: struct, simulation parameters
%--------------------------------------------------------------------------
% Simulate firm characteristics
real_param.alpha = sim_param.theta_0*data.zj_alpha+sqrt(sim_param.Lambda_0)*randn(sim_param.J,1);
real_param.beta = sim_param.theta_1*data.zj_beta+sqrt(sim_param.Lambda_1)*randn(sim_param.J,1);
real_param.tau = sim_param.psi*data.zj_tau+sqrt(sim_param.delta)*randn(sim_param.J,1);
real_param.v = exp(real_param.tau);
% Simulate returns data
data.f_t = sim_param.mu_f + sqrt(sim_param.Omega_f)*randn(sim_param.T,1);
data.regressor = [ones(sim_param.T,1),data.f_t];
data.y_jt = ones(sim_param.T,1)*real_param.alpha'+data.f_t*real_param.beta'+...
    repmat(sqrt(real_param.v'),sim_param.T,1).*randn(sim_param.T,sim_param.J);
end