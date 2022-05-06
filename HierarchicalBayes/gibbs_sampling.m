function est = gibbs_sampling(initial,data,sim_param)
%--------------------------------------------------------------------------
% Estimate the model parameters with Gibbs sampleing
% initial: struct, preliminary estimates from least squares
% data: struct, market data of stock returns and factor returns
% sim_param: struct, simulation parameters
%--------------------------------------------------------------------------
param_names = fieldnames(initial);
for i = 1:length(param_names)
    est.(param_names{i}) = zeros(size(initial.(param_names{i})));
end
param_samples = initial;
for g = 1:sim_param.G
    param_samples = gibbs_iteration(param_samples,data,sim_param);
    if g > sim_param.B
        for i = 1:length(param_names)
            est.(param_names{i}) = est.(param_names{i})+...
                param_samples.(param_names{i})/(sim_param.G-sim_param.B);
        end
    end
end
end