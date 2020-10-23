function Y=blockdiag_unfold_cell(X,r_n,n1,n2,n3)

r=cumsum(r_n);
r=[0 r];

[m,n]=size(X);
Y=cell(1,n3);



for i=0:1:n3-1
    Y{i+1}=X(n1*i+1:n1*(i+1),r(i+1)+1:r(i+2));
end

end