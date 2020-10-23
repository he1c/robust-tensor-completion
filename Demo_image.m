%%----------------demo video-----------------

clear

addpath(genpath(pwd))

I=imread('data\dance.jpg');

I=double(I)./255;

[n1,n2,n3]=size(I);

%% algorithm settings

option=[];

option.debug    = 0;
option.maxitr   = 500;
option.stopc    = 1e-5;

% HQ method
option.sigmamin = 0.15;
option.qtmin    = 0.25;
option.yita     = 2;

% rank
option.rank     = 100;

%% system settings

p  = 0.5;
v1 = 0.001;
v2 = 1;
c = 0.2;

%% initial

omega = find(rand(n1*n2*n3,1)<p);
Mask = zeros(n1,n2,n3);
Mask(omega) = 1;

G = zeros(n1,n2,n3);
for i=1:1:n3
    G(:,:,i)=noisemix(n1,n2,c,v1,v2,'gaussian');
end
I_n=I+G;

MissM=Mask.*I_n;

%% HQ-HQ_TCASD
tic
option.stopc  = 1e-7;
option.stopc2 = 1e-5;
option.lambda = 0.2;
F_HQ_TCASD = HQ_TCASD(MissM,Mask,I,option);
err = F_HQ_TCASD-I;
disp(['HQ_TCASD : ' num2str(10*log10(n1*n2*n3/norm(err(:))^2)) ' , time : ' num2str(toc)])

%% HQ-TCTF
tic
option.stopc   = 1e-5;
option.method  = 'image';
F_HQ_TCTF = HQ_TCTF(MissM,Mask,I,option);
err = F_HQ_TCTF-I;
disp(['HQ_TCTF : ' num2str(10*log10(n1*n2*n3/norm(err(:))^2)) ' , time : ' num2str(toc)])
    


