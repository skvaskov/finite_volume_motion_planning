function [meanQ,stdQ] = compute_mean_and_stddev(Q,m,grid_lower_bounds,grid_upper_bounds)
d = length(m);


h = (grid_lower_bounds-grid_upper_bounds)./m;


xgrid = get_grid_points(m,grid_lower_bounds,grid_upper_bounds);

meanQ = NaN(length(Q),d);
stdQ = NaN(length(Q),d);

for i = 1:d
    
    other_idxs = (1:d)~=i;
    
    meanQ(:,i) = sum(((xgrid(:,i)+h(i)/2).^2/2-(xgrid(:,i)-h(i)/2).^2/2).*cat(2,Q{:})*prod(h(other_idxs)))';
    
    stdQ(:,i) = sqrt(sum(((xgrid(:,i)+h(i)/2).^3/3-(xgrid(:,i)-h(i)/2).^3/3).*cat(2,Q{:})*prod(h(other_idxs)))'-meanQ(:,i).^2);
end


end

% function [meanQ,stdQ] = compute_mean_and_stddev(Q,d,m)
% 
% x = linspace(0,1,m+1);
% 
% h = 1/m;
% 
% sub_i = ind2subnd(m*ones(1,d),1:m^d);
% 
% xgrid = x(sub_i)+h/2;
% 
% meanQ = NaN(length(Q),d);
% stdQ = NaN(length(Q),d);
% 
% for i = 1:d
% meanQ(:,i) = sum(repmat(((xgrid(:,i)+h/2).^2/2-(xgrid(:,i)-h/2).^2/2),[1 length(Q)]).*cat(2,Q{:})*h^(d-1))';
% stdQ(:,i) = sqrt(sum(repmat(((xgrid(:,i)+h/2).^3/3-(xgrid(:,i)-h/2).^3/3),[1 length(Q)]).*cat(2,Q{:})*h^(d-1))'-meanQ(:,i).^2);
% end
% 
% 
% end