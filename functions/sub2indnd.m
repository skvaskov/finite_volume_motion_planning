%submatrix is Nxd

function ind_i = sub2indnd(SZ,sub_matrix)

d = size(sub_matrix,2);

nonzero_indexj_s = mat2cell(sub_matrix',ones(1,d));

        
ind_i = sub2ind(SZ,nonzero_indexj_s{:});
        
ind_i = ind_i';
end