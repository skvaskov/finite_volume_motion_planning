%inputs Q: 1xNtimestep cell array. Each cell contains prod(m)x1 vector with grid values at each cell
        %dproj vector containing the dimensions you want to project onto
        %m: dx1 number of grid cells in each dimension
        %grid_lower(upper)_bounds: dx1 lower (upper) limit of grid in each dimenstion
        %assume the grid is rectangular with uniform spacing
        
%outputs: Qproj: 1xNtimestep cell array. Each cell contains prod(m(dproj))x1 vector with grid values at each cell

%written by Sean Vaskov on May 11th, 2020

function Qproj = project_onto_dimension(Q,dproj,m,grid_lower_bounds,grid_upper_bounds)

d = length(m);

other_idxs = ~ismember((1:d),dproj);

h = (grid_upper_bounds-grid_lower_bounds)./m;


Qproj = cell(size(Q));


szproj = m(dproj);

sub_i = ind2subnd(m,1:prod(m));
idx = sub2indnd(szproj,sub_i(:,dproj));

ProjMat = sparse(idx,1:prod(m),prod(h(other_idxs)),prod(m(dproj)),prod(m));

for n = 1:length(Q)
    
    
    Qproj{n} = ProjMat*Q{n};
    
    
end


end

