function Qproj = project_onto_dimension(Q,dproj,m,grid_lower_bounds,grid_upper_bounds)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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

