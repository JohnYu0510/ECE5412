%% ===================================================================
% Simulation Study on the Price Limits Model with Gibbs Sampling
% Jinyuan Yu (jy478) 2022.04
%% ===================================================================
close all; clear; clc; rng(5412);
%% -------------------------------------------------------------------
% Set Global Parameters
% -------------------------------------------------------------------------
% Basic parameters
sim_param.T = 1000; % Number of days for estimation
% Market characteristics
sim_param.mu_f = 0.005; % Mean of the market index return
sim_param.Omega_f = 0.04^2/20; % Variance of the market index return
% Firm characteristics
sim_param.beta = [0;1]; % Factor loadings
sim_param.sigma2 = 0.08^2; % Idiosyncratic risk
% Gibbs samplings parameters
sim_param.G = 2000; % Number of Gibbs iteration
sim_param.B = 500; % Number of discarded samples
% Prior Belief
prior.b = [0.01;1.5];prior.B = [0.001,0.0001;0.0001,0.001];
prior.delta = 0.005^2*sim_param.T; prior.v = 0.5;
prior.b = [0;0];prior.B = [0,0;0,0];
prior.v = 1; prior.delta = 0.1;
% Price Limits
sim_param.lu = 0.1;
sim_param.ld = -0.1;

%% -------------------------------------------------------------------
% Simulate Market Data
% -------------------------------------------------------------------------
% Initialize datasets
data.lu = sim_param.lu; data.ld = sim_param.ld; data.T = sim_param.T;
% Simulate factor returns
data.ft = nan(sim_param.T,2);
data.ft(:,1) = 1;
data.ft(:,2) = sim_param.mu_f+sqrt(sim_param.Omega_f)*randn(sim_param.T,1);
% Simulate true stock returns
real.r = data.ft*sim_param.beta+sqrt(sim_param.sigma2)*randn(sim_param.T,1);
% Compute censored returns
E = 0; data.z = nan(sim_param.T,1);
for t = 1:sim_param.T
    if real.r(t,1)+E>sim_param.lu
        data.z(t,1) = sim_param.lu;
        E = (real.r(t,1)+E-sim_param.lu);
    elseif real.r(t,1)+E<sim_param.ld
        data.z(t,1) = sim_param.ld;
        E = (real.r(t,1)+E-sim_param.ld);
    else
        data.z(t,1) = real.r(t,1)+E;
        E = 0;
    end
end
data.real_r = real.r;
% Cluster the observed returns
data.group_idx = cluster_returns(data.z,data.lu,data.ld);

%% -------------------------------------------------------------------
% Gibbs sampling
% -------------------------------------------------------------------------
% Initialize the iteration
ls.beta = inv(data.ft'*data.ft)*data.ft'*data.z;
ls.sigma2 = (data.z-data.ft*ls.beta)'*(data.z-data.ft*ls.beta)/data.T;
input = ls; input.r = data.real_r;%initial_r(ls.beta,ls.sigma2,data);
gibbs.beta = nan(2,sim_param.G);
gibbs.sigma2 = nan(1,sim_param.G);
gibbs.r = nan(data.T,sim_param.G);

% Start iteration
for g = 1:sim_param.G
    output = gibbs_iteration(input,data,prior);
    gibbs.beta(:,g) = output.beta;
    gibbs.sigma2(1,g) = output.sigma2;
    gibbs.r(:,g) = output.r;
    input = output;
    if mod(g,100) == 0
        fprintf('%.1f\n',g)
    end
end

%% -------------------------------------------------------------------
% Print results
% -------------------------------------------------------------------------
display(mean(gibbs.beta(:,sim_param.B+1:sim_param.G),2))
display(mean(gibbs.sigma2(:,sim_param.B+1:sim_param.G),2))