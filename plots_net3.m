%script to compare MC & Analy

which_set=1; %change to see either pulse (1) or sinusoidal input (2)

switch which_set %assuming same time vectors (tmv)
    case 1 %fast pulse input
        load dQmc_n3Tw
        load dAn_n3Tw
        load dMC_n3Tw
    
    case 2 %sinusoidal input
        load dQmc_n3Tw_sin
        load dAn_n3Tw_sin
        load dMC_n3Tw_sin
end


%specify all (or subset) of Nc for mean/var, and Nc*(Nc-1)/2 for covars
ind_mv=(1:Nc)';

ind_cv=[];
for j=1:Nc
    for k=j+1:Nc
        ind_cv=[ind_cv; k j];
    end
end
subSetI=ones(Nc*(Nc-1)*.5,1); %rand(Nc*(Nc-1)*.5,1)<0.3;
Jin=subSetI.*ind_cv(:,1);
Kin=subSetI.*ind_cv(:,2);

cc=parula(length(ind_mv)+1);
cc1=parula(length(ind_cv)+1);

figure
hold on
for j=1:length(ind_mv)
    plot(tmv,mnX_Mqmc(j,:),'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,mnX_M(j,:),'k--','LineWidth',3)
    plot(tmv,mnX_Ma(j,:),'color',cc(j,:))
end
set(gca,'FontSize',26)
xlabel('Time (a.u.)')
ylabel('Mean activity')

figure
hold on
for j=1:length(ind_mv)
    plot(tmv,mnF_Mqmc(j,:),'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,mnF_M(j,:),'k--','LineWidth',3)
    plot(tmv,mnF_Ma(j,:),'color',cc(j,:))
end
set(gca,'FontSize',26)
ylabel('Mean firing')

figure
hold on
for j=1:length(ind_mv)
    vTmp=squeeze(covX_M(j,j,:)); %cov X
    vTmpa=squeeze(covX_Ma(j,j,:)); %cov X
    vTmpqm=squeeze(covX_Mqmc(j,j,:)); %cov X
    
    plot(tmv,vTmpqm,'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,vTmp,'k--','LineWidth',3)
    plot(tmv,vTmpa,'color',cc(j,:))
end
set(gca,'FontSize',26)
ylabel('Var activity')

figure
hold on
for j=1:length(ind_mv)
    vTmp=squeeze(covF_M(j,j,:)); %cov F
    vTmpa=squeeze(covF_Ma(j,j,:)); %cov F
    vTmpqm=squeeze(covF_Mqmc(j,j,:)); %cov F

    plot(tmv,vTmpqm,'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,vTmp,'k--','LineWidth',3)
    plot(tmv,vTmpa,'color',cc(j,:))
end
set(gca,'FontSize',26)
ylabel('Var firing')

figure
hold on
for j=1:length(ind_cv)
    
    cvTmp=squeeze(covX_M(Jin(j),Kin(j),:)); %cov X
    cvTmpa=squeeze(covX_Ma(Jin(j),Kin(j),:)); %cov X
    cvTmpqm=squeeze(covX_Mqmc(Jin(j),Kin(j),:)); %cov X
    
    plot(tmv,cvTmpqm,'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,cvTmp,'k--','LineWidth',3)
    plot(tmv,cvTmpa,'color',cc1(j,:))
    
end
set(gca,'FontSize',26)
ylabel('Cov activity')

figure
hold on
for j=1:length(ind_cv)
    cvTmp=squeeze(covF_M(Jin(j),Kin(j),:)); %cov X
    cvTmpa=squeeze(covF_Ma(Jin(j),Kin(j),:)); %cov X
    cvTmpqm=squeeze(covF_Mqmc(Jin(j),Kin(j),:)); %cov X
    
    plot(tmv,cvTmpqm,'LineWidth',1.5,'color',[.5 .5 .5])
    plot(tmv,cvTmp,'k--','LineWidth',3)
    plot(tmv,cvTmpa,'color',cc1(j,:))
end
set(gca,'FontSize',26)
ylabel('Cov firing')


