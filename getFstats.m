function [mn_Fa,cov_Fa]=getFstats(mnXa,covXa,rv_vec,sp_vec)
%[mn_Fa,cov_Fa]=getFstats(mnXa,covXa,rv_vec,sp_vec);
%gets firing rate stats from mean activity
%OUTPUTS: 
% mn_Fa = Nc x 1 vector of mean firing rates
% cov_Fa = Nc x Nc covariance matrix of firing rate activity

Nc=length(mnXa);

%need discret norm pdfs to do analytic calcs
d_nrm=0.01;
norm_dom=(-3:d_nrm:3)';
len_nrm=length(norm_dom);
normDom_mat=repmat(norm_dom,1,Nc); %matrx len_nrm x Nc, use in F_mn/vr calcs
norm01_pdf=normpdf(norm_dom,0,1); %stand. normal
nrm_varF=repmat(norm_dom,len_nrm,1); %[x1(1),x1(2),...]
nrm_varS=reshape(repmat(norm_dom',len_nrm,1),len_nrm*len_nrm,1); %[x2(1),x2(1),...]    

norm2D_sing=zeros(len_nrm*len_nrm,1);
rv_Mat=repmat(rv_vec',len_nrm,1);
sp_Mat=repmat(sp_vec',len_nrm,1);

corr_Xa=zeros(Nc,Nc); %entire correl matrix of Xj (ones on diagonal)
%OUTPUTS; iterative solver, start at uncoupled values
mn_Fa=zeros(Nc,1);    %mean of F(Xj)
cov_Fa=zeros(Nc,Nc);  %entire cov matrix of Fj (includes var)



var_t=diag(covXa);
varX_matr=var_t*var_t'; %used to divide to get correl matrix
%calc the corr(Xj,Xk)
corr_Xa=covXa./sqrt(varX_matr);

%check if correl matrix is positive definite (assuming semi-def is unlikely)
[L,p_chol]=chol(corr_Xa);
if(p_chol>0)
    disp(['Corr/cov is not PSD.']);
    %return; %exit function; corr structure is invalid
end

%calc firing stats, assuming (X_j,Y_k) are multivariate Normal
Finst_1D=0.5*(1+tanh((normDom_mat*diag(sqrt(var_t))+repmat(mnXa',len_nrm,1)-rv_Mat)./sp_Mat));
mn_Fa=norm01_pdf'*Finst_1D*d_nrm; %mean firing rate
mn_Fa=mn_Fa'; %make it Nc x 1

%set the variance of the firing rate
cov_Fa=diag( (norm01_pdf'*( Finst_1D.^2 )*d_nrm)'-mn_Fa.^2 );

%Recalc numeric 2D-Gauss with FINAL cov/corr, calc all .5Nc(Nc-1) of them
%not saving 2D-Gauss; finish calculation here, cov_Fa off-diagonal
for colInd=2:Nc
    for rowInd=1:(colInd-1)
        if(corr_Xa(rowInd,colInd)>=1)
            sigTmp=[1 1-1e-5; 1-1e-5  1];
        elseif(corr_Xa(rowInd,colInd)<=-1)
            sigTmp=[1 -1+1e-5; -1+1e-5  1];
        else
            sigTmp=[1 corr_Xa(rowInd,colInd); corr_Xa(rowInd,colInd)  1];
        end
        [L,p_chol]=chol(sigTmp);

        norm2D_sing=mvnpdf([nrm_varF nrm_varS],[0 0],sigTmp); %0 mean, unit var, cov=corr_Xa
        
        cov_Fa(rowInd,colInd)=...
            sum(repmat(Finst_1D(:,rowInd),len_nrm,1).*reshape(repmat(Finst_1D(:,colInd)',len_nrm,1),len_nrm*len_nrm,1)...
            .*norm2D_sing)*d_nrm*d_nrm - mn_Fa(rowInd)*mn_Fa(colInd);
        
        cov_Fa(colInd,rowInd)=cov_Fa(rowInd,colInd); %make symmetric
    end
end


