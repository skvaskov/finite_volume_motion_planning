%input: f: function Z \times K to R^d that gives the trajectory dynamics, d
%          is the state dimension
%         where z is the state vector and k are parameters
%         n_k: number of parameters
%         h: cell widths in each dimension
%output: fbar: 1xd cell array. each cell contains a function handle the
%        inputs to this function is a 1xd array containing the left cell boundary
%        and a 1xn_k vector with the trajectory parameters k. The output of
%        this function is the average velocity at the cell interface

%created by Sean Vaskov on May 17th, 2020

function fbar = compute_average_velocity_functions(f,n_k,h)
    d = length(h);
    
    fbar = cell(1,d);
    
    %create symbolic variables
    z = sym('z',[1 d]);
    xlow = sym('x',[1 d]);
    k = sym('k',[1 n_k]);
    
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

