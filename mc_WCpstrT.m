function [cov_F,mn_F,cov_X,mn_X]=mc_WCpstrT(Nc,Pstruct,t,mu_fn,sig_fn)
%[cov_F,mn_F,cov_X,mn_X]=mc_WCpstrT(Nc,Pstruct,t,mu_t,sig_fn);
%Simulate full network of Nc (input) WC, coupled/corrNoise, etc
%sig_vec, tau_vec are all Nc x 1 column vectors!
%F is sigmoidal; 1/2*(1+tanh((Inp-rv_vec)./sp_vec))
%Gm is an NcxNc coupling matrix, with no autaptic (diag(Gm)=0)
%CinMat is an NcxNc correlation matrix; ones on diag and PSD!
% last 2 inputs: t, mu_t determine the time-varying input to the
% cells; t is Lt x 1, mu_t is an Nc X Lt matrix
%--- should check: (sig_vec,tau_vec,rv_vec,sp_vec) are all Nc x 1
%---    CinMat is PSD, Gm is Nc x Nc ---
% OUTPUT: cov_F,cov_X are NcXNcXLt, matrices, mn_F,mn_X are NcXLt matrices

rng('shuffle') %seed random number generator

% Other parameters
% Extract parameters
tau_vec     = Pstruct.tau_vec;
rv_vec      = Pstruct.rv_vec;
sp_vec      = Pstruct.sp_vec;
Gm          = Pstruct.Gm;
CinMat      = Pstruct.CinMat;

%%!! assuming sig_vec constant!! otherwise, must change R_cor (noise) at
%%all steps!
sig_vec = sig_fn(0); 

if(length(sig_vec)~=Nc || length(tau_vec)~=Nc || length(rv_vec)~=Nc || length(sp_vec)~=Nc ... 
    || size(Gm,1)~=Nc || size(Gm,2)~=Nc || size(CinMat,1)~=Nc)
    disp('Check params!');
    return;
end

dt=t(2)-t(1);
Lt=length(t);
sq_dt=1/sqrt(dt);


N_relz=1000000;
if(N_relz*Nc > 1e8) %limit to how large matrices can be
    disp('Size & #Realz too large!');
    return;
end
    
tauInv=1./tau_vec;
itauMat=repmat(tauInv',N_relz,1);
covM_xi=CinMat.*(sig_vec*sig_vec'); %pure cov matrix
R_cor=chol(covM_xi); %R_cor'*randn(Nc,1) gives correct noise term (RHS)

rvMat=repmat(rv_vec',N_relz,1);
spMat=repmat(sp_vec',N_relz,1);

%Start sims
muT   = mu_fn(0);
xS=randn(N_relz,Nc)*diag(sig_vec./sqrt(2*tau_vec))+repmat(muT',N_relz,1); %# of cells
%do a t=5 run to get correct IC starting point
for k=1:round(5/dt)
    xi=randn(N_relz,Nc);
    FxM = 0.5*(1+tanh((xS-rvMat)./spMat));
    muT   = mu_fn(t(1));  %Use current time
    %eqns for cells 
    xS=xS+dt*itauMat.*( -xS+repmat(muT',N_relz,1)+sq_dt*(xi*R_cor)+FxM*(Gm') );
end

%OUTPUTS
mn_X=zeros(Nc,Lt);    %mean of Xj
cov_X=zeros(Nc,Nc,Lt);  %entire cov matrix of Xj (includes var)
mn_F=zeros(Nc,Lt);    %mean of F(Xj)
cov_F=zeros(Nc,Nc,Lt);  %entire cov matrix of Fj (includes var)

%in main loop, since N_relz x Nc, use (L*xi')'=xi*L', where L'=R_cor
% similarly for (Gm*FxMat')'=FxMat*Gm'

for k=1:Lt
    xi=randn(N_relz,Nc);
    FxM = 0.5*(1+tanh((xS-rvMat)./spMat));
    
    muT   = mu_fn(t(k));  %Use current time
    
    %eqns for cells 
    xS=xS+dt*itauMat.*( -xS+repmat(muT',N_relz,1)+sq_dt*(xi*R_cor)+FxM*(Gm') );
    
    %keep running sum of stats/dens
    mn_X(:,k)=mean(xS)';
    mn_F(:,k)=mean(FxM)';
    cov_X(:,:,k)=xS'*xS./N_relz - mn_X(:,k)*mn_X(:,k)';
    cov_F(:,:,k)=FxM'*FxM./N_relz - mn_F(:,k)*mn_F(:,k)';
    
end




