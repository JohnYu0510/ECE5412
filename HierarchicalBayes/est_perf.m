function performance = est_perf(est_param,real_param)
%--------------------------------------------------------------------------
% Evaluate the performance of estimation
% est_param: struct, estimated parameters
% real_param, struct, real parameters
%--------------------------------------------------------------------------
performance = table;
performance.alpha_mae = mean(abs(est_param.alpha-real_param.alpha));
performance.beta_mae = mean(abs(est_param.beta-real_param.beta));
performance.variance_mae = mean(abs(est_param.v-real_param.v));
end