function F_o=HQ_TCASD(MissM,Mask,I,option)

%% HQ-TCASD
sigmamin = option.sigmamin;
stopc = option.stopc;
stopc2 = option.stopc2;
qtmin = option.qtmin;
debug = option.debug;
maxitr = option.maxitr;
yita = option.yita;
lambda = option.lambda;

[n1,n2,n3]=size(MissM);

r=option.rank;

Mhat=fft(MissM,[],3);

for i=1:1:n3
    X{i}=fft(0.001*rand(n1,r),[],3);
    Y{i}=fft(0.001*rand(r,n2),[],3);
    M{i}=Mhat(:,:,i);
end

F_o=zeros(n1,n2,n3);

n_MissM=norm(MissM(:))^2;

err=100;

%clear Mhat

diff_count=0;

sigma=100;

for kk=1:1:maxitr
    
    %% kernetl width selection
       
    J=Mask.*(MissM-F_o);
    AAA=J(:);
    AAA(AAA==0)=[];
    
    if sigma~=sigmamin   
        sigma=max(max(abs(quantile(AAA,qtmin)),abs(quantile(AAA,1-qtmin)))*yita,sigmamin);   
    end   

    MaskMCC=Mask.*exp(-J.^2/sigma^2/2);
         
    %% update X
    
    Res=cellfun(@(B,C,D)(B-C*D),M,X,Y,'UniformOutput',false);
    
    dX_M=fft(MaskMCC.*ifft(cat(3,Res{:}),[],3),[],3);

    dX=unfold_cell(dX_M);

    dX_blk=cellfun(@(A,B)A*B',dX,Y,'UniformOutput',false);

    dX_blk2=cellfun(@(A,B)A*B,dX_blk,Y,'UniformOutput',false);

    tx_2 = sqrt(MaskMCC).*ifft(cat(3,dX_blk2{:}),[],3);

    temp=cellfun(@(A)norm(A,'fro'),dX_blk,'UniformOutput',false);
        
    ss=cat(3,temp{:});
       
    tx_1=sum(ss.^2)/n3;
    
    tx= tx_1 / norm(tx_2(:))^2;
       
    for j=1:1:n3
        X{j}=X{j}+tx*dX_blk{j};
    end   
    
    %clear tx_2 dX_blk dX_blk2 dX_M dX 
     
    %% update Y

    Res=cellfun(@(B,C,D)(B-C*D),M,X,Y,'UniformOutput',false);

    dY_M=fft(MaskMCC.*ifft(cat(3,Res{:}),[],3),[],3)/n3;

    dY=unfold_cell(dY_M);

    dY=cellfun(@(A,B)A'*B,X,dY,'UniformOutput',false);
    
    dY_blk2=cellfun(@(A,B)A*B,X,dY,'UniformOutput',false);

    ty1_1 =cat(3,dY{:});

    ty1_2 = sqrt(MaskMCC).*ifft(cat(3,dY_blk2{:}),[],3);

    ty1= norm(ty1_1(:))^2 / norm(ty1_2(:))^2;
   
    dY2=cellfun(@(A,B)pinv(A'*A)*B*n3,X,dY,'UniformOutput',false);

    dY12=cellfun(@(A,B)conj(A).*B,dY,dY2,'UniformOutput',false);

    dY_blk2=cellfun(@(A,B)A*B,X,dY2,'UniformOutput',false);

    ty2_1 =cat(3,dY12{:});

    ty2_2 = sqrt(MaskMCC).*ifft(cat(3,dY_blk2{:}),[],3);

    ty2 = abs(sum(ty2_1(:))) / norm(ty2_2(:))^2;

    for j=1:1:n3
        Y{j}=Y{j}+(1-lambda)*ty1*dY{j}+lambda*ty2*dY2{j};
    end
    
    %clear ty1_2 ty2_2 dY_M dY dY2 dY_blk2

    %% Output and converge check

    FF=cellfun(@(A,B)A*B,X,Y,'UniformOutput',false);

    C_pre=sqrt(MaskMCC).*(MissM-F_o); 

    F_o=real(ifft(cat(3,FF{:}),[],3));

    C=sqrt(MaskMCC).*(MissM-F_o);
    
    err_pre=err;
    
    err=(norm(C(:))^2-norm(C_pre(:))^2)/n_MissM;
    
    diff=err-err_pre;  % add stop ceriterion by checking the second order of error (optional)
    
    if abs(diff)<stopc||abs(err)<stopc2
        diff_count=diff_count+1;
        if diff_count>1  % avoid early terminate
            break;
        end
    end
   
    if debug
        err1=I-F_o;
        disp(norm(err1(:)));
    end

 end

end
