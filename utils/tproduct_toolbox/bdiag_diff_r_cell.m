function Y=bdiag_diff_r_cell(X)

[~,n3]=size(X);

Y=[];

for i=1:1:n3
    Y = blkdiag(Y,X{i});
end


end