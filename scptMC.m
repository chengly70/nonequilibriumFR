%script to run corresponding monte carlo AFTER running scptAn_[].m

%get params from Analy run
load('dAn_n3Tw.mat','Nc','Pstruct','tmv','mu_fn','sig_fn')
%load('dAn_n3Tw_sin.mat','Nc','Pstruct','tmv','mu_fn','sig_fn') %for sinus input
%--duplicate previous line to LOAD Nc=50, etc. etc.)--

Lt=length(tmv);

% -- outputs to save --
mnX_M=zeros(Nc,Lt);
covX_M=zeros(Nc,Nc,Lt);
mnF_M=zeros(Nc,Lt);
covF_M=zeros(Nc,Nc,Lt);

tic
    [covF_M,mnF_M,covX_M,mnX_M]=mc_WCpstrT(Nc,Pstruct,tmv,mu_fn,sig_fn); 
toc

save dMC_n3Tw mnX_M covX_M mnF_M covF_M tmv
%save dMC_n3Tw_sin mnX_M covX_M mnF_M covF_M tmv %for sinusoidal input
%--duplicate previous line to SAVE Nc=50, etc. etc.)--