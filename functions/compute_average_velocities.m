%input: k: 1xn_k array containing the trajectory parameters k
%       fbar: 1xd cell array containing function handles for the average
%       velocity in each dimension at the cell interfaces, output of
%       compute_average_velocity_functions
%       m: 1xd array with the number of cells in each dimension
%       grid_lower(upper)_bounds: 1xd array with the lower(upper) bounds of
%       the grid. assume rectangular and unfirom spacing
%
%output: ubar: dx1 cell array. each cell contains a vectory with the
%        function fbar evaluated at the cell interfaces in each dimension

%Created by Sean Vaskov on May 17th, 2020

function ubar = compute_average_velocities(k,fbar,m,grid_lower_bounds,grid_upper_bounds)

d = length(m);

h = (grid_upper_bounds-grid_lower_bounds)./m;

ubar = cell(d,1);


for j = 1:d
    
    other_idxs = (1:d)~=j;
        
    nj = (m(j)+1)*prod(m(other_idxs));
    
    szd = m;
    szd(j) = m(j)+1;
    
    sub_i = ind2subnd(szd,1:nj);
    lowerlims = (sub_i-1).*repmat(h,[nj,1])+grid_lower_bounds;

    
    ubar{j} = fbar{j}(lowerlims,k);


end 

end

