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

