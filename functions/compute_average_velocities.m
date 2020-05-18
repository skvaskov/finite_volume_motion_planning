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

