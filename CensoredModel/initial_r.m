function r = initial_r(beta,sigma2,data)
%--------------------------------------------------------------------------
% Generate initial r for iteration
%--------------------------------------------------------------------------
r = nan(size(data.z)); 
for t = 1:size(r,1)
    if data.z(t)==data.lu
        pd = makedist('normal','mu',data.ft(t,:)*beta,'sigma',sqrt(sigma2));
        r(t,1) = random(truncate(pd,data.lu,inf),1,1);
    elseif data.z(t)==data.ld
        pd = makedist('normal','mu',data.ft(t,:)*beta,'sigma',sqrt(sigma2));
        r(t,1) = random(truncate(pd,-inf,data.ld),1,1);
    else
        r(t,1) = data.z(t,1);
    end
end