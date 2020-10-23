function F_o=TCASD(MissM,Mask,I,option)

stopc=option.stopc;
debug=option.debug;
maxitr=option.maxitr;
lambda=option.lambda;

[n1,n2,n3]=size(MissM);

r=option.rank;

Mhat=fft(MissM,[],3);

for i=1:1:n3
    XX{i}=abs(0.001*rand(n1,r)); % 0.001 for image and video, 1 for synthetic data
    YY{i}=abs(0.001*rand(r,n2)); % 0.001 for image and video, 1 for synthetic data
    X{i}=fft(XX{i},[],3);
    Y{i}=fft(YY{i},[],3);
    M{i}=Mhat(:,:,i);
end

F_o=zeros(n1,n2,n3);

n_MissM=norm(MissM(:))^2;

for kk=1:1:maxitr
       
    %% update X

    Res=cellfun(@(B,C,D)(B-C*D),M,X,Y,'UniformOutput',false);

    Res_M = cat(3,Res{:});

    dX_M=fft(Mask.*ifft(Res_M,[],3),[],3);

    dX=unfold_cell(dX_M);

    dX_blk=cellfun(@(A,B)A*B',dX,Y,'UniformOutput',false);

    dX_blk2=cellfun(@(A,B)A*B,dX_blk,Y,'UniformOutput',false);

    tx_2 = Mask.*ifft(cat(3,dX_blk2{:}),[],3);

    temp=cellfun(@(A)norm(A,'fro'),dX_blk,'UniformOutput',false);
    ss=cat(3,temp{:});
    tx_1=sum(ss.^2)/n3;

    tx= tx_1 / norm(tx_2(:))^2;

    for j=1:1:n3
        X{j}=X{j}+tx*dX_blk{j};
    end
    
    %% update Y

    Res=cellfun(@(B,C,D)(B-C*D),M,X,Y,'UniformOutput',false);

    Res_M = cat(3,Res{:});

    dY_M=fft(Mask.*ifft(Res_M,[],3),[],3)/n3;

    dY=unfold_cell(dY_M);

    dY=cellfun(@(A,B)A'*B,X,dY,'UniformOutput',false);
    
    dY_blk2=cellfun(@(A,B)A*B,X,dY,'UniformOutput',false);

    ty1_1 =cat(3,dY{:});

    ty1_2 = Mask.*ifft(cat(3,dY_blk2{:}),[],3);

    ty1= norm(ty1_1(:))^2 / norm(ty1_2(:))^2;
   
    dY2=cellfun(@(A,B)pinv(A'*A)*B*n3,X,dY,'UniformOutput',false);

    dY12=cellfun(@(A,B)conj(A).*B,dY,dY2,'UniformOutput',false);

    dY_blk2=cellfun(@(A,B)A*B,X,dY2,'UniformOutput',false);

    ty2_1 =cat(3,dY12{:});

    ty2_2 = Mask.*ifft(cat(3,dY_blk2{:}),[],3);

    ty2 = abs(sum(ty2_1(:))) / norm(ty2_2(:))^2;

    for j=1:1:n3
        Y{j}=Y{j}+(1-lambda)*ty1*dY{j}+lambda*ty2*dY2{j};
    end
    
    %% Output and converge check

    FF=cellfun(@(A,B)A*B,X,Y,'UniformOutput',false);

    C_pre=Mask.*(MissM-F_o); 

    F_o=real(ifft(cat(3,FF{:}),[],3));

    C=Mask.*(MissM-F_o);
    
    err=(norm(C(:))^2-norm(C_pre(:))^2)/n_MissM;
    
    if abs(err)<stopc
        break;
    end
    
    if debug
        err1=I-F_o;
        disp(norm(err1(:)))
    end

end

end
