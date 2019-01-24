%script to compare MC & Analy
%now with sin inputs

cc1=hsv(4);

for whSet=1:4 %loop over all 4 Gm strengths
    
switch whSet %assuming same time vectors (tmv)
    case 1 %weakest coupling
        load(['dAn_n50sin_a.mat'],'tmv','Nc','mnX_Ma','mnF_Ma','covF_Ma','covX_Ma')
        load(['dMC_n50sin_a.mat'],'mnX_M','mnF_M','covF_M','covX_M')
    case 2
        load(['dAn_n50sin_b.mat'],'tmv','Nc','mnX_Ma','mnF_Ma','covF_Ma','covX_Ma')
        load(['dMC_n50sin_b.mat'],'mnX_M','mnF_M','covF_M','covX_M')
    case 3
        load(['dAn_n50sin_c.mat'],'tmv','Nc','mnX_Ma','mnF_Ma','covF_Ma','covX_Ma')
        load(['dMC_n50sin_c.mat'],'mnX_M','mnF_M','covF_M','covX_M')
    case 4 %strongest coupling
        load(['dAn_n50sin_d.mat'],'tmv','Nc','mnX_Ma','mnF_Ma','covF_Ma','covX_Ma')
        load(['dMC_n50sin_d.mat'],'mnX_M','mnF_M','covF_M','covX_M')
end
%L1-norm of error
Err_mnX=abs(mnX_Ma-mnX_M);
Err_mnF=abs(mnF_Ma-mnF_M);
Err_cvX=abs(covX_Ma-covX_M);
Err_cvF=abs(covF_Ma-covF_M);

figure(142)
hold on
plot(tmv,mean(Err_mnX),'color',cc1(whSet,:),'LineWidth',2) 
%errorbar(tmv,mean(Err_mnX),std(Err_mnX),'color',cc1(whSet,:)) 
plot(tmv,mean(Err_mnF),'-.','color',cc1(whSet,:)) 
set(gca,'FontSize',22)
xlabel('Time (a.u.)')
ylabel('Error in mean')


figure(144)
hold on
evTmp=zeros(size(tmv));
for j=1:Nc
    evTmp=evTmp+squeeze(Err_cvX(j,j,:)); %err var of X
end
evTmp=evTmp./Nc; %average L1-Error
plot(tmv,evTmp,'color',cc1(whSet,:),'LineWidth',2) 

figure(145)
hold on
evTmp=zeros(size(tmv));
for j=1:Nc
    evTmp=evTmp+squeeze(Err_cvF(j,j,:)); %err var of X
end
evTmp=evTmp./Nc; %average L1-Error
plot(tmv,evTmp,'color',cc1(whSet,:),'LineWidth',2) 


figure(144)
hold on
ecvTmp=zeros(size(tmv));
for j=1:Nc
    for k=j+1:Nc
    ecvTmp=ecvTmp+squeeze(Err_cvX(j,k,:)); %cov X
    
    end
end
ecvTmp=ecvTmp./(Nc*(Nc-1)*.5);
plot(tmv,ecvTmp,'color',cc1(whSet,:))
set(gca,'FontSize',22)
xlabel('Time (a.u.)')
ylabel('Error covar activity')


figure(145)
hold on
ecvTmp=zeros(size(tmv));
for j=1:Nc
    for k=j+1:Nc
    ecvTmp=ecvTmp+squeeze(Err_cvF(j,k,:)); %cov X
    end
end
ecvTmp=ecvTmp./(Nc*(Nc-1)*.5);
plot(tmv,ecvTmp,'color',cc1(whSet,:))
set(gca,'FontSize',22)
xlabel('Time (a.u.)')
ylabel('Error covar firing')

end

