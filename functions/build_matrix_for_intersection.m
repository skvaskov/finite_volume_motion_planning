function M = build_matrix_for_intersection(m,d,ftprint)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
szq = m*ones(1,d);

sub_i = ind2subnd(szq,1:m^d);


Mr = cell(m^d,1);
Ir = cell(m^d,1);

for i = 1:m^d
    Lr = (sub_i(:,1)-sub_i(i,1) <= ftprint(1)) | (sub_i(:,2)-sub_i(i,2) <= ftprint(2));
    
    Mr{i} = find(Lr);
    Ir{i} = i*ones(nnz(Lr),1);
end


end

