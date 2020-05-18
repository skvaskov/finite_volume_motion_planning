%inputs x:Nxd values of states
        %m: widthr of grid cells in each dimension
        %grid_lower_bounds: lower limit of grid in each dimenstion
        %assume the grid is rectangular with uniform spacing
        
%outputs: sub_out: Nxd indexs of grid cells in each dimension      

function sub_out = x_val_to_sub(x,h,grid_lower_bounds)

N = size(x,1);

sub_out = floor((x-grid_lower_bounds)./repmat(h,[N 1]))+1;

end
