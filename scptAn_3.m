%script to test method on Nc=3

Nc=3; %larger network
tau_vec=1+0.1*randn(Nc,1);
mu_vec=rand(Nc,1) - 0.5; %[-.5,.5]
sig_vec=rand(Nc,1) + 1;

rv_vec=0.1*randn(Nc,1);
sp_vec=0.35*rand(Nc,1) + 0.05;
Cin_sub=0.8*randn(Nc,Nc);
Cin_sub=Cin_sub'*Cin_sub; %make symme
scl_dia=diag(1./sqrt(diag(Cin_sub))); %diag matrix to put 1's diag CinMat
CinMat=scl_dia*Cin_sub*scl_dia;  %should be positive definite, check:
R=chol(CinMat);

%len_vr=4; %# times vary params
vr=1*.1;
Gm=zeros(Nc,Nc); %vary coupling matrix, randomly chosen
%change Gm or ..?
gm_tmp=vr*randn(Nc,Nc);
%Gm(:,:,j)=gm_tmp-diag(diag(gm_tmp));
Gm=gm_tmp;

dt=0.01;
tmv=(0:dt:8)'; %arbitrary time vector
Lt=length(tmv);

% -- outputs to save --
mnX_Ma=zeros(Nc,Lt);
covX_Ma=zeros(Nc,Nc,Lt);
mnF_Ma=zeros(Nc,Lt);
covF_Ma=zeros(Nc,Nc,Lt);


% Make this a function call instead
Ton=2;
tauR=0.2;
tauD=0.5;
Io=0.5;
Imx=1.75;
mu_fn=@(t) Io+(Imx-Io)*(t>Ton).*(exp(-(t-Ton)./tauD)-exp(-(t-Ton)./tauR))./(tauD-tauR)*ones(Nc,1);
% For sigma: make (possibly) time-dependent
sig_fn = @(t) sig_vec;

% Set up parameter structure
Pstruct.tau_vec     = tau_vec;
Pstruct.rv_vec  = rv_vec;
Pstruct.sp_vec  = sp_vec;
Pstruct.Gm      = Gm ;
Pstruct.CinMat  = CinMat;

% Set up options
options = odeset('RelTol',1e-5,'AbsTol',1e-8);
%set IC
mu_vBG   = mu_fn(0);
sig_vec  = sig_fn(0);
[convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa]=iter_method(Nc,mu_vBG,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);
mnX_Ma(:,1)=mn_Xa;
covX_Ma(:,:,1)=cov_Xa;
mnF_Ma(:,1)=mn_Fa;
covF_Ma(:,:,1)=cov_Fa;
% Set up IC for call to ODE45
initC    = [mn_Xa; reshape(cov_Xa+mn_Xa*mn_Xa',Nc^2,1)];


tic

[tout,yout] = ode45(@rhs_momentsFP,tmv,initC,options,Pstruct,mu_fn,sig_fn);


% get the firing stats here
for ind_Sn=2:Lt
    %removed your abs() in covX_Ma that was causing a sign-error
    mnX_Ma(:,ind_Sn)  = yout(ind_Sn,1:Nc);
    covX_Ma(:,:,ind_Sn)=reshape(yout(ind_Sn,Nc+1:end),Nc,Nc)-mnX_Ma(:,ind_Sn)*mnX_Ma(:,ind_Sn)';
    
    [mn_Fa,cov_Fa]=getFstats(mnX_Ma(:,ind_Sn),covX_Ma(:,:,ind_Sn),rv_vec,sp_vec);
    mnF_Ma(:,ind_Sn)=mn_Fa;
    covF_Ma(:,:,ind_Sn)=cov_Fa;
end


toc

save dAn_n3Tw covF_Ma covX_Ma mnX_Ma mnF_Ma Nc Pstruct vr tmv mu_fn sig_fn
