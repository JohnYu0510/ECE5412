function group_idx = cluster_returns(z,lu,ld)
%--------------------------------------------------------------------------
% Create groups of consecutive limit days and a subsequent non-limit day
% z: censored observed returns data
%--------------------------------------------------------------------------
loc_pos = 1; group_idx = nan(size(z,1),2); num_rows = 0;
for t = 1:size(z,1)
    if ~((z(t,1)==lu)||(z(t,1)==ld))
        group_idx(num_rows+1,:) = [loc_pos,t];
        loc_pos = t+1;
        num_rows = num_rows+1;
    end
group_idx = group_idx(1:num_rows,:);
end