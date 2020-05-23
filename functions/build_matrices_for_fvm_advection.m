%input: m: 1xd array containing the number of cells in each axis
%       ubar: 1xd cell array. each cell contains a (m(d)+1)prod(m(~=d)) x d
%       vector containing the velocity at the i-1/2 cell interface
        
%outputs: [UP,DQI,DQi,DF] matrices that are used int fvm_advection function

%written by Sean Vaskov on May 16th, 2020

function [UP,DQI,DQi,DF] = build_matrices_for_fvm_advection(ubar,m)

d = length(ubar);

UP = cell(d,1);
DF = cell(d,1);
DQI = cell(d,1);
DQi = cell(d,1);

for j = 1:d
    
    
    nj = (m(j)+1)*prod(m((1:d)~=j));
    
    UP_idx  = NaN(nj,1);

    DQI_pos_idx = NaN(nj,1);
    DQI_neg_idx = NaN(nj,1);
    
    szd = m;
    szd(j) = m(j)+1;
    
    zzd = zeros(1,d);
    zzd(j) = 1;
    
    sub_ih = ind2subnd(szd,(1:nj));
    sub_ih_minus_1 = sub_ih - repmat(zzd,[nj 1]);
    sub_ih_minus_2 = sub_ih - 2*repmat(zzd,[nj 1]);
    sub_ih_plus_1  = sub_ih + repmat(zzd,[nj 1]);
    
    %clip at 1,m
    sub_ih(sub_ih(:,j) < 1, j) = 1;
    sub_ih(sub_ih(:,j) > m(j), j) = m(j);
    
    sub_ih_minus_1(sub_ih_minus_1(:,j) < 1, j) = 1;
    sub_ih_minus_1(sub_ih_minus_1(:,j) > m(j), j) = m(j);
    
    sub_ih_minus_2(sub_ih_minus_2(:,j) < 1, j) = 1;
    sub_ih_minus_2(sub_ih_minus_2(:,j) > m(j), j) = m(j);
    
    sub_ih_plus_1(sub_ih_plus_1(:,j) < 1, j) = 1;
    sub_ih_plus_1(sub_ih_plus_1(:,j) > m(j), j) = m(j);
    
    ind_i = sub2indnd(m,sub_ih);
    ind_i_minus_1 = sub2indnd(m,sub_ih_minus_1);
    ind_i_minus_2 = sub2indnd(m,sub_ih_minus_2);
    ind_i_plus_1 = sub2indnd(m,sub_ih_plus_1);
    
    Lp = ubar{j} >= 0;
    Ln = ~Lp;
    
    DQI_pos_idx(Lp) = ind_i_minus_1(Lp);
    DQI_neg_idx(Lp) = ind_i_minus_2(Lp);
    UP_idx(Lp) = ind_i_minus_1(Lp);
    
    DQI_pos_idx(Ln) = ind_i_plus_1(Ln);
    DQI_neg_idx(Ln) = ind_i(Ln);
    UP_idx(Ln) = ind_i(Ln);
    

    UP{j} = sparse(1:nj,UP_idx,ubar{j},nj,prod(m));
    DQI{j} = sparse([1:nj,1:nj],[DQI_pos_idx;DQI_neg_idx],[ones(nj,1);-ones(nj,1)],nj,prod(m));
    DQi{j} = sparse([1:nj,1:nj],[ind_i;ind_i_minus_1],[ones(nj,1);-ones(nj,1)],nj,prod(m));
    
    sub_i = ind2subnd(m,1:prod(m));
    sub_p = sub_i;
    sub_p(:,j) = sub_p(:,j)+1;
    
    ind_i = sub2indnd(szd,sub_i);
    ind_p = sub2indnd(szd,sub_p);
   
    DF{j} = sparse([1:prod(m),1:prod(m)],[ind_p;ind_i],[ones(prod(m),1);-ones(prod(m),1)],prod(m),nj);

end


end

