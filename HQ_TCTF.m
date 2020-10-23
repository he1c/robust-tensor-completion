function TC = HQ_TCTF( MissM,Mask,I,option )

%% produce data
Omega=find(Mask(:)==1);
data=MissM(Omega);
known=Omega;
r_n=option.rank;

[n1,n2,n3]=size(MissM);
Nway=[n1,n2,n3];

%% our method
if strcmp(option.method,'fixrank')
    option.rank_adj = 0*ones(1,n3);
    option.rank_min = r_n*ones(1,n3);
    option.EstCoreNway = r_n*ones(1,n3);
elseif strcmp(option.method,'adjust')
    option.rank_adj = -1*ones(1,n3);
    option.rank_min = 5*ones(1,n3);
    option.EstCoreNway = round(r_n*ones(1,n3));
elseif strcmp(option.method,'image')
    option.rank_adj = -1*ones(1,n3);
    option.rank_min = [100,25*ones(1,n3-1)];
    option.EstCoreNway = round(option.rank*ones(1,n3));
elseif strcmp(option.method,'video')
    option.rank_adj = -1*ones(1,n3);
    option.rank_min = [80,60*ones(1,n3-1)];
    option.EstCoreNway = round(option.rank*ones(1,n3));
else
    disp('error : unkown type of data')
    return;
end

TC = HQ_TCTF_solver(MissM,I,Mask,option);

end

