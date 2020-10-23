function Y=T_LSX(T,X,A)

%% construct B
[m,n,k]=size(T);
[~,r,~]=size(X);
T_t=fft(T,[],3)*k;
X_t=fft(X,[],3);
X_l=bdiag(X_t);

parfor j=1:n

    b=reshape(T_t(:,j,:),[],1);
    
    AX=A{j}*X_l;   
       
    AA=AX'*AX;
  
    y_l=AA\(AX'*b);
    
    Y_l(:,j,:)=reshape(y_l,r,k);
       
end 

    
Y=ifft(Y_l,[],3);

Y=real(Y);

end