function fbar = compute_average_velocity_functions(f,d_k,h)
    d = length(h);
    
    fbar = cell(1,d);
    
    %create symbolic variables
    z = sym('z',[1 d]);
    xlow = sym('x',[1 d]);
    k = sym('k',[1 d_k]);
    
    fsym = f(z,k);
    
    %create function handle for each dimension
    for i = 1:d
        
        other_idxs = find((1:d)~=i);
        
        ftemp = fsym(i);
        
        %take average over other dimensions
        for j = other_idxs
            
            ftemp = int(ftemp,z(j),xlow(j),xlow(j)+h(j))/h(j);
            
        end
        
        ftemp = subs(ftemp,z(i),xlow(i));
        
        fbar{i} = matlabFunction(ftemp,'vars',{xlow,k});
        
        
    end
    
    
    
end

