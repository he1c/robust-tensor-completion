function TC = HQ_TCTF_solver(T,I,Mask,opts)

%% Data preprocessing and initialization
[n1,n2,n3]=size(T);
rank_adj = opts.rank_adj;
rank_min = opts.rank_min;
sigmamin = opts.sigmamin;
maxitr= opts.maxitr;
stopc = opts.stopc;
qtmin = opts.qtmin;
debug = opts.debug;
yita = opts.yita;
coreNway = opts.EstCoreNway;

%% initialization X and Y
for i=1:1:n3
    XX{i}=abs(0.001*rand(n1,coreNway(i)));
    YY{i}=abs(0.001*rand(coreNway(i),n2));
    X{i}=fft(XX{i},[],3);
    Y{i}=fft(YY{i},[],3);
end

%% compute the initialization residual
C=zeros(n1,n2,n3);
for n = 1:n3
    C(:,:,n)=X{n}*Y{n};
end

rho=0.95;
mu=1.01;

TC=ifft(C,[],3);

normT=norm(T(:));

sigma=100;

TC_pre=10000;

for kk = 1:maxitr
    
  
    %% update kernel width (modified)

    J=Mask.*(T-TC);
    AAA=J(:);
    AAA(AAA==0)=[];
    
    if sigma~=sigmamin  
        sigma=max(max(abs(quantile(AAA,qtmin)),abs(quantile(AAA,1-qtmin)))*yita,sigmamin);    
    end

    MCC=Mask.*exp(-(Mask.*(T-TC)).^2/sigma^2/2);

    TC=TC+MCC./(1+MCC).*Mask.*(T-TC); 
    
    C=fft(TC,[],3);

    %% update (X,Y)  (same as TCTF)
    for n = 1:n3
        X{n}=C(:,:,n)*Y{n}'*pinv(Y{n}*Y{n}');
        Xsq{n} = X{n}'*X{n};
        Y{n} = pinv(Xsq{n})*X{n}'*C(:,:,n);
        C(:,:,n)=X{n}*Y{n};
    end
    
    %% adjust the rank of (X,Y) (same as TCTF)
    if rank_adj(n) == -1 && rho<1
        max_k=max(coreNway);
        sum_k=sum(coreNway);
        sigmas=zeros(max_k*n3,1);
        for i=1:n3
            s = svd(Xsq{i});
            sigmas((i-1)*max_k+1:(i-1)*max_k+length(s))=s;
        end
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);
        [dmx,imx] = max(drops);
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);
        if rel_drp>10
            thold=rho*sum(dR);
            iidx=0;ss=0;
            len=length(dR);
            for i=1:len
                ss=ss+dR(i);
                if(ss>thold)
                    iidx=i;
                    break;
                end
            end
            if(iidx>sum(rank_min(n)))
                idx=floor((id(iidx+1:sum_k)-1)/max_k);
                for n=1:n3
                    num=length(find(idx==n-1));
                    if(num>0)
                        if coreNway(n)-num>rank_min(n)
                            coreNway(n) = coreNway(n)-num;
                        else
                            coreNway(n) = rank_min(n);
                        end
                        [Qx,Rx] = qr(X{n},0);
                        [Qy,Ry] = qr(Y{n}',0);
                        [U,S,V] = svd(Rx*Ry');
                        sigv = diag(S);
                        X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
                        Y{n} = (Qy*V(:,1:coreNway(n)))';
                        C(:,:,n)=X{n}*Y{n};
                    end
                end
            end
            rho=rho*mu;
        end
    end
    
    
    %% judge whether converges (modified)
    
    MaskMCCsq=sqrt(MCC);
    
    CC_pre=MaskMCCsq.*(T-TC_pre);
    
    TC=real(ifft(C,[],3));
   
    CC=MaskMCCsq.*(T-TC); 
    
    TC_pre=TC;

    err=(norm(CC(:))^2-norm(CC_pre(:))^2)/normT^2;
    
    if kk>10&&abs(err)<stopc
        break;
    end
    
    if debug
        err1=I-TC;
        disp(norm(err1(:)))
    end
    
end

end