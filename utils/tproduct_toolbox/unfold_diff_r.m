function Y=unfold_diff_r(X)

[~,n3]=size(X);

Y=[];

for i=1:1:n3
    Y = [Y;X{i}];
end


end