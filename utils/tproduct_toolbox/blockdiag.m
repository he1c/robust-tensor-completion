function Y=blockdiag(X,r_n,n1,n2,n3)

[m,n]=size(X);
Y=zeros(m,n);

r=cumsum(r_n);
r=[0 r];

for i=0:1:n3-1
    Y(n1*i+1:n1*(i+1),r(i+1)+1:r(i+2))=X(n1*i+1:n1*(i+1),r(i+1)+1:r(i+2));
end

end