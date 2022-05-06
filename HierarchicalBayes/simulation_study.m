%% ===================================================================
% Simulation Study on the Hierarchical Bayesian Factor Model
% Jinyuan Yu (jy478) 2022.04
%% ===================================================================
close all; clear; clc; rng(5412);
%% -------------------------------------------------------------------
% Set Global Parameters
% -------------------------------------------------------------------------
% Basic parameters
sim_param.J = 30; % Number of stocks
sim_param.T = 24; % Number of months for estimation
sim_param.M = 100; % Number of simulation repliacations
% Market characteristics
sim_param.mu_f = 0.01; % Mean of the market index return
sim_param.Omega_f = 0.04^2; % Variance of the market index return
% Firm characteristics
data.zj_alpha = ones(sim_param.J,1);
data.zj_beta = ones(sim_param.J,1);
data.zj_tau = ones(sim_param.J,1);
% Hierarchical Bayesian model hyperparameters
sim_param.theta_0 = 0; sim_param.Lambda_0 = 0.007^2; 
sim_param.theta_1 = 1; sim_param.Lambda_1 = 0.0025^2; 
sim_param.delta = 4*(log(0.052^2+0.087^2)-log(0.087^2));
sim_param.psi = 4*log(0.097)-log(0.054^2+0.087^2);
% Gibbs samplings parameters
sim_param.G = 2000; % Number of Gibbs iteration
sim_param.B = 1000; % Number of discarded samples

%% -------------------------------------------------------------------
% Simulate Market Data
% -------------------------------------------------------------------------
% Initialize the estimation accuracy performance
ls_perf = table; hb_perf = table;
% Run the estimation
for m = 1:sim_param.M
    % Simulate data
    [real_param,data] = simulate_data(data,sim_param);
    % Initailize the estimation with least squares
    initial = preliminary_estimate(data,sim_param);
    % Estimate with Gibbs sampling
    hb_est = gibbs_sampling(initial,data,sim_param);
    % Performance
    ls_perf = [ls_perf;est_perf(initial,real_param)];
    hb_perf = [hb_perf;est_perf(hb_est,real_param)];
    fprintf('%4.f\n',m)
end
%% -------------------------------------------------------------------
% Print results
% -------------------------------------------------------------------------
fprintf('alpha_MAE,beta_MAE,variance_MAE\n')
fprintf('ls: %.4f,%.4f,%.4f\n',mean(ls_perf.alpha_mae),mean(ls_perf.beta_mae),mean(ls_perf.variance_mae))
fprintf('hb: %.4f,%.4f,%.4f\n',mean(hb_perf.alpha_mae),mean(hb_perf.beta_mae),mean(hb_perf.variance_mae))