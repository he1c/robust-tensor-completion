%%----------------demo-----------------

clear

addpath(genpath(pwd))

%% experiment settings
n1 = 200;
n2 = 200;
n3 = 20;
r = 10; % tubal rank
p = 0.5; % observation percentage

%% GMM noise setting
v1 = 0.01;
v2 = 4;
c = 0.1;

%% algorithm settings

option=[];

option.debug    = 0;
option.maxitr   = 500;
option.stopc    = 1e-10; 

% HQ method
option.sigmamin = 0.3;
option.qtmin    = 0.25;
option.yita     = 6;

% rank
option.rank     = r;


%% Synthetic data
Tx=randn(n1,r,n3);
Ty=randn(r,n2,n3);
I = tprod(Tx,Ty); % size: n1*n2*n3

G = zeros(n1,n2,n3);
for i=1:1:n3
    G(:,:,i)=noisemix(n1,n2,c,v1,v2,'gaussian');
end
X=I+G;

omega = find(rand(n1*n2*n3,1)<p);
MissM = zeros(n1,n2,n3);
Mask = MissM;
Mask(omega) = 1;
MissM(omega) = X(omega);

%% TCASD
option.stopc  = 1e-9;
option.lambda = 0;
tic
F_TCASD=TCASD(MissM,Mask,I,option);
err=F_TCASD-I;
disp(['TCASD : ' num2str(norm(err(:))/norm(I(:))) ' , time : ' num2str(toc)])

%% HQ-TCASD (lambda=0)
tic
option.lambda = 0;
option.stopc  = 1e-11; 
option.stopc2 = 1e-9;
F_HQ_TCASD = HQ_TCASD(MissM,Mask,I,option);
err = F_HQ_TCASD-I;
disp(['HQ-TCASD-0 : ' num2str(norm(err(:))/norm(I(:))) ' , time : ' num2str(toc)])

%% HQ-TCASD (lambda=1)
tic
option.lambda = 1;
F_HQ_TCASD_S = HQ_TCASD(MissM,Mask,I,option);
err = F_HQ_TCASD_S-I;
disp(['HQ-TCASD-l : ' num2str(norm(err(:))/norm(I(:))) ' , time : ' num2str(toc)])

%% HQ-TCTF
tic
option.method   = 'fixrank';
option.stopc = 1e-9; 
F_HQ_TCTF = HQ_TCTF(MissM,Mask,I,option);
err = F_HQ_TCTF-I;
disp(['HQ-TCTF : ' num2str(norm(err(:))/norm(I(:))) ' , time : ' num2str(toc)])




