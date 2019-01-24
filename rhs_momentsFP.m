function [rhs_FP] = rhs_momentsFP(t,momentsFP,Pstruct,mu_fn,sig_fn)
% RHS for ODE for 1st and 2nd moments of time-dependent Fokker-Planck
%  
%  INPUTS: 
%     t:            time
%     momentsFP:    (Nc + Nc^2) x 1 array  : 1st and 2nd moments (not
%       centered)
%     Pstruct:      structure w/ all parameters (except for TD mu, sigma
%     mu_fn:       function pointer returning: Nc x 1 array, taking a
%               single real argument (time)
%     sig_fn:       function pointer returning: Nc x 1 array, taking a
%               single real argument (time)
%
% OUTPUTS:
%     rhs_FP:   RHS for ODE 

% Extract, reshape moments
Nc  =  (sqrt(4*length(momentsFP)+1)-1)/2;   % Must assume correct size!

mu0   = momentsFP(1:Nc);
Ejk0  = reshape(momentsFP(Nc+1:end),Nc,Nc);

Var0  = diag(Ejk0)-mu0.*mu0;

% Extract parameters
tau_vec     = Pstruct.tau_vec;
rv_vec      = Pstruct.rv_vec;
sp_vec      = Pstruct.sp_vec;
Gm          = Pstruct.Gm;
CinMat      = Pstruct.CinMat;

mu_vec    = mu_fn(t);
sig_vec   = sig_fn(t);

% Call RHS functions

% Note: rhs_mn takes Var0
Mn_rhs      = rhs_mn(mu0,mu_vec,Var0,tau_vec,rv_vec,sp_vec,Gm);
Ejk_rhs     = rhs_cov1(Ejk0,mu0,mu_vec,Var0,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);

if (~isreal(Mn_rhs))
    warning(sprintf('Mn_rhs contains imaginary part: t = %g',t));
end

% Package for return
rhs_FP = [Mn_rhs; reshape(Ejk_rhs,Nc^2,1)];

end

