function output = gibbs_iteration(input,data,prior)
%--------------------------------------------------------------------------
% One step of Gibbs iterations to update the parameter estimators
% input: struct, parameters simulated from the last iteration
% data: struct, market data of stock returns and factor returns
% prior: struct, to characterize the prior distribution of parameters
%--------------------------------------------------------------------------
% Generate a posterior sample for beta
Omega_hat = inv((prior.B+(data.ft'*data.ft)/input.sigma2));
beta_hat = Omega_hat*(prior.B*prior.b+(data.ft'*input.r)/input.sigma2);
U = chol(Omega_hat);
output.beta = beta_hat+U'*randn(size(beta_hat,1),1);
% Generate a posterior sample for sigma2
term1 = (prior.v+data.T)/2; 
term2 = prior.delta+(input.r-data.ft*output.beta)'*(input.r-data.ft*output.beta)/2;
%pd = makedist('Gamma','a',term1,'b',1/term2);
output.sigma2 = 1/gamrnd(term1,1/term2);
% Generate posterior samples for true returns
output.r = nan(data.T,1);
for i = 1:size(data.group_idx,1)
    Ti = data.group_idx(i,2) - data.group_idx(i,1) + 1;
    zi = data.z(data.group_idx(i,1):data.group_idx(i,2));
    ri = input.r(data.group_idx(i,1):data.group_idx(i,2));
    fi = data.ft(data.group_idx(i,1):data.group_idx(i,2),:);
    if Ti == 1
        output.r(data.group_idx(i,2),1) = data.z(data.group_idx(i,2),1);
        continue
    end
    i_lb = nan(Ti-1,1); i_ub = nan(Ti-1,1);
    for j = 1:(Ti-1)
        j_lb_k = -inf(Ti-j,1); j_ub_k = inf(Ti-j,1);
        for k = j:(Ti-1)
            if (k==1) Esum = 0; else Esum = sum(ri(1:(k-1))-zi(1:(k-1))); end
            if (zi(k,1) == data.lu)
                j_lb_k(k-j+1,1) = data.lu - ri(k,1) - Esum + ri(j);
                j_ub_k(k-j+1,1) = nan;
            else
                j_ub_k(k-j+1,1) = data.ld- ri(k,1) - Esum + ri(j);
                j_lb_k(k-j+1,1) = nan;
            end
        end
        j_lb = max(j_lb_k); j_ub = min(j_ub_k);
        if isnan(j_lb) j_lb = -inf; end
        if isnan(j_ub) j_ub = inf; end
        pd = makedist('normal','mu',fi(j,:)*output.beta,'sigma',sqrt(output.sigma2));
        ri(j,1) = random(truncate(pd,j_lb,j_ub),1,1);
        i_lb(j) = j_lb; i_ub(j) = j_ub;
    end
    ri(Ti,1) = sum(zi) - sum(ri(1:(Ti-1),:));
    output.r(data.group_idx(i,1):data.group_idx(i,2)) = ri;
end
end