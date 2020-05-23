%input: SZ: 1xd array contining the matrix size
%    sub_i: 1xd array containing the position in each dimension for the
%        cell indexed by i

%        
%output: i: double containing the index (position in column vector)

%Created by Sean Vaskov on March 4th, 2020

function ind_i = sub2indnd(SZ,sub_matrix)

d = size(sub_matrix,2);

nonzero_indexj_s = mat2cell(sub_matrix',ones(1,d));

        
ind_i = sub2ind(SZ,nonzero_indexj_s{:});
        
ind_i = ind_i';
end