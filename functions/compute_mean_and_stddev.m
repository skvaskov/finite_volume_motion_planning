%input: Q: 1xNtimestep cell array. Each cell contains a prod(m)x1 array
%          containing the cell values 
%       m: 1xd vector with the number of cells in each dimension
%       grid_lower(upper)_bounds: 1xd lower (upper) bound of state space in each dimension
        
%outputs: meanQ: Ntimestep x d array containing the average value in each
%                dimension at each timestep
%         StdQ: Ntimestep x d array containing the standard deviation in each
%                dimension at each timestep

%written by Sean Vaskov on May 15th, 2020

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
