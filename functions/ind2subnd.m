%input: SZ: 1xd array contining the matrix size
%        i: double containing the index (position in column vector)
%output: sub_i: 1xd array containing the position in each dimension for the
%        cell indexed by i

%Created by Sean Vaskov on March 4th, 2020

function sub_i = ind2subnd(SZ,i)

%transpose if mutliple entries and given vertically
if size(i,1) ~= 1 && size(i,2) == 1
    i = i';
end

d = size(SZ,2);
        sub_i = cell(d,1);
        [sub_i{:}] = ind2sub(SZ,i);
        sub_i = cell2mat(sub_i)';

end

