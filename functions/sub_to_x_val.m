%inputs sub_in:Nxd indices of cell along each axis
        %h: width of grid cells in each dimension
        %grid_lower_bounds: lower limit of grid in each dimenstion
        %assume the grid is rectangular with uniform spacing
        
%outputs: x_out: Nxd grid point corresponding to cell center     

function x_out = sub_to_x_val(sub_in,h,grid_lower_bounds)

N = size(sub_in,1);

x_out = (sub_in-1/2).*repmat(h,[N 1])+grid_lower_bounds;



end
