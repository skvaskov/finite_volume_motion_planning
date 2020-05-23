%input: m: 1xd vector with the number of cells in each dimension
%       grid_lower(upper)_bounds: 1xd lower (upper) bound of state space in each dimension
        
%outputs: xgrid: prod(m) x d double containing the center coordinate of
%each cell

%written by Sean Vaskov on May 16th, 2020

function xgrid = get_grid_points(m,grid_lower_bounds,grid_upper_bounds)
d = length(m);

N = prod(m);

h = (grid_upper_bounds-grid_lower_bounds)./m;

sub_i = ind2subnd(m,1:N);

xgrid = NaN(N,d);

for i = 1:d
    x = linspace(grid_lower_bounds(i),grid_upper_bounds(i),m(i)+1);
    xgrid(:,i) = x(sub_i(:,i))'+h(i)/2;
end


end

