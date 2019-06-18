Norm_D_90 = dlmread('500x8__Normal_D_90.txt','\t');
Norm_Ds_90 = dlmread('500x8__Normal_Ds_90.txt','\t');
Norm_L_90 = dlmread('500x8__Normal_LinSpace_90.txt','\t');
SSA_D_90 = dlmread('500x8__SSA_D_90.txt','\t');
SSA_Ds_90 = dlmread('500x8__SSA_Ds_90.txt','\t');
SSA_L_90 = dlmread('500x8__SSA_LinSpace_90.txt','\t');

Norm_D_60 = dlmread('500x8__Normal_D_60.txt','\t');
Norm_Ds_60 = dlmread('500x8__Normal_Ds_60.txt','\t');
Norm_L_60 = dlmread('500x8__Normal_LinSpace_60.txt','\t');
SSA_D_60 = dlmread('500x8__SSA_D_60.txt','\t');
SSA_Ds_60 = dlmread('500x8__SSA_Ds_60.txt','\t');
SSA_L_60 = dlmread('500x8__SSA_LinSpace_60.txt','\t');

Norm_D_120 = dlmread('500x8__Normal_D_120.txt','\t');
Norm_Ds_120 = dlmread('500x8__Normal_Ds_120.txt','\t');
Norm_L_120 = dlmread('500x8__Normal_LinSpace_120.txt','\t');
SSA_D_120 = dlmread('500x8__SSA_D_120.txt','\t');
SSA_Ds_120 = dlmread('500x8__SSA_Ds_120.txt','\t');
SSA_L_120 = dlmread('500x8__SSA_LinSpace_120.txt','\t');

theta_t = [0.5,3,9,3];

ndm_90=[];
ndsm_90=[];
nlm_90=[];
sdm_90=[];
sdsm_90=[];
slm_90=[];

for i=1:4
    ndm_90=[ndm_90 mean(Norm_D_90(:,i))];
    ndsm_90=[ndsm_90 mean(Norm_Ds_90(:,i))];
    nlm_90=[nlm_90 mean(Norm_L_90(:,i))];
    sdm_90=[sdm_90 mean(SSA_D_90(:,i))];
    sdsm_90=[sdsm_90 mean(SSA_Ds_90(:,i))];
    slm_90=[slm_90 mean(SSA_L_90(:,i))];
end

ndm_60=[];
ndsm_60=[];
nlm_60=[];
sdm_60=[];
sdsm_60=[];
slm_60=[];

for i=1:4
    ndm_60=[ndm_60 mean(Norm_D_60(:,i))];
    ndsm_60=[ndsm_60 mean(Norm_Ds_60(:,i))];
    nlm_60=[nlm_60 mean(Norm_L_60(:,i))];
    sdm_60=[sdm_60 mean(SSA_D_60(:,i))];
    sdsm_60=[sdsm_60 mean(SSA_Ds_60(:,i))];
    slm_60=[slm_60 mean(SSA_L_60(:,i))];
end

ndm_120=[];
ndsm_120=[];
nlm_120=[];
sdm_120=[];
sdsm_120=[];
slm_120=[];

for i=1:4
    ndm_120=[ndm_120 mean(Norm_D_120(:,i))];
    ndsm_120=[ndsm_120 mean(Norm_Ds_120(:,i))];
    nlm_120=[nlm_120 mean(Norm_L_120(:,i))];
    sdm_120=[sdm_120 mean(SSA_D_120(:,i))];
    sdsm_120=[sdsm_120 mean(SSA_Ds_120(:,i))];
    slm_120=[slm_120 mean(SSA_L_120(:,i))];
end

%%
biases_90 = [norm(ndm_90-theta_t);norm(ndsm_90-theta_t);norm(nlm_90-theta_t);...
    norm(sdm_90-theta_t);norm(sdsm_90-theta_t);norm(slm_90-theta_t)];
covariances_90=[det(cov(Norm_D_90));det(cov(Norm_Ds_90));det(cov(Norm_L_90));...
    det(cov(SSA_D_90));det(cov(SSA_Ds_90));det(cov(SSA_L_90))];
sMSEd_90 = [0 0 0 0];
for i=1:size(SSA_D_90,1)
    sMSEd_90=sMSEd_90+(theta_t - SSA_D_90(i,:)).*(theta_t - SSA_D_90(i,:))/size(SSA_D_90,1);
end
sMSEds_90 = [0 0 0 0];
for i=1:size(SSA_Ds_90,1)
    sMSEds_90=sMSEds_90+(theta_t - SSA_Ds_90(i,:)).*(theta_t - SSA_Ds_90(i,:))/size(SSA_Ds_90,1);
end
sMSEl_90 = [0 0 0 0];
for i=1:size(SSA_L_90,1)
    sMSEl_90=sMSEl_90+(theta_t - SSA_L_90(i,:)).*(theta_t - SSA_L_90(i,:))/size(SSA_L_90,1);
end
nMSEd_90 = [0 0 0 0];
for i=1:size(Norm_D_90,1)
    nMSEd_90=nMSEd_90+(theta_t - Norm_D_90(i,:)).*(theta_t - Norm_D_90(i,:))/size(Norm_D_90,1);
end
nMSEds_90 = [0 0 0 0];
for i=1:size(Norm_D_90,1)
    nMSEds=nMSEds_90+(theta_t - Norm_Ds_90(i,:)).*(theta_t - Norm_Ds_90(i,:))/size(Norm_Ds_90,1);
end
nMSEl_90 = [0 0 0 0];
for i=1:size(Norm_L_90,1)
    nMSEl_90=nMSEl_90+(theta_t - Norm_L_90(i,:)).*(theta_t - Norm_L_90(i,:))/size(Norm_L_90,1);
end

mse_90 = [norm(sMSEd_90);norm(sMSEds_90);norm(sMSEl_90);norm(nMSEd_90);norm(nMSEds_90);norm(nMSEl_90)];
%%
biases_60 = [norm(ndm_60-theta_t);norm(ndsm_60-theta_t);norm(nlm_60-theta_t);...
    norm(sdm_60-theta_t);norm(sdsm_60-theta_t);norm(slm_60-theta_t)];
covariances_60=[det(cov(Norm_D_60));det(cov(Norm_Ds_60));det(cov(Norm_L_60));...
    det(cov(SSA_D_60));det(cov(SSA_Ds_60));det(cov(SSA_L_60))];
sMSEd_60 = [0 0 0 0];
for i=1:size(SSA_D_60,1)
    sMSEd_60=sMSEd_60+(theta_t - SSA_D_60(i,:)).*(theta_t - SSA_D_60(i,:))/size(SSA_D_60,1);
end
sMSEds_60 = [0 0 0 0];
for i=1:size(SSA_Ds_60,1)
    sMSEds_60=sMSEds_60+(theta_t - SSA_Ds_60(i,:)).*(theta_t - SSA_Ds_60(i,:))/size(SSA_Ds_60,1);
end
sMSEl_60 = [0 0 0 0];
for i=1:size(SSA_L_60,1)
    sMSEl_60=sMSEl_60+(theta_t - SSA_L_60(i,:)).*(theta_t - SSA_L_60(i,:))/size(SSA_L_60,1);
end
nMSEd_60 = [0 0 0 0];
for i=1:size(Norm_D_60,1)
    nMSEd_60=nMSEd_60+(theta_t - Norm_D_60(i,:)).*(theta_t - Norm_D_60(i,:))/size(Norm_D_60,1);
end
nMSEds_60 = [0 0 0 0];
for i=1:size(Norm_D_60,1)
    nMSEds=nMSEds_60+(theta_t - Norm_Ds_60(i,:)).*(theta_t - Norm_Ds_60(i,:))/size(Norm_Ds_60,1);
end
nMSEl_60 = [0 0 0 0];
for i=1:size(Norm_L_60,1)
    nMSEl_60=nMSEl_60+(theta_t - Norm_L_60(i,:)).*(theta_t - Norm_L_60(i,:))/size(Norm_L_60,1);
end

mse_60 = [norm(sMSEd_60);norm(sMSEds_60);norm(sMSEl_60);norm(nMSEd_60);norm(nMSEds_60);norm(nMSEl_60)];
%%
biases_120 = [norm(ndm_120-theta_t);norm(ndsm_120-theta_t);norm(nlm_120-theta_t);...
    norm(sdm_120-theta_t);norm(sdsm_120-theta_t);norm(slm_120-theta_t)];
covariances_120=[det(cov(Norm_D_120));det(cov(Norm_Ds_120));det(cov(Norm_L_120));...
    det(cov(SSA_D_120));det(cov(SSA_Ds_120));det(cov(SSA_L_120))];
sMSEd_120 = [0 0 0 0];
for i=1:size(SSA_D_120,1)
    sMSEd_120=sMSEd_120+(theta_t - SSA_D_120(i,:)).*(theta_t - SSA_D_120(i,:))/size(SSA_D_120,1);
end
sMSEds_120 = [0 0 0 0];
for i=1:size(SSA_Ds_120,1)
    sMSEds_120=sMSEds_120+(theta_t - SSA_Ds_120(i,:)).*(theta_t - SSA_Ds_120(i,:))/size(SSA_Ds_120,1);
end
sMSEl_120 = [0 0 0 0];
for i=1:size(SSA_L_120,1)
    sMSEl_120=sMSEl_120+(theta_t - SSA_L_120(i,:)).*(theta_t - SSA_L_120(i,:))/size(SSA_L_120,1);
end
nMSEd_120 = [0 0 0 0];
for i=1:size(Norm_D_120,1)
    nMSEd_120=nMSEd_120+(theta_t - Norm_D_120(i,:)).*(theta_t - Norm_D_120(i,:))/size(Norm_D_120,1);
end
nMSEds_120 = [0 0 0 0];
for i=1:size(Norm_D_120,1)
    nMSEds=nMSEds_120+(theta_t - Norm_Ds_120(i,:)).*(theta_t - Norm_Ds_120(i,:))/size(Norm_Ds_120,1);
end
nMSEl_120 = [0 0 0 0];
for i=1:size(Norm_L_120,1)
    nMSEl_120=nMSEl_120+(theta_t - Norm_L_120(i,:)).*(theta_t - Norm_L_120(i,:))/size(Norm_L_120,1);
end

mse_120 = [norm(sMSEd_120);norm(sMSEds_120);norm(sMSEl_120);norm(nMSEd_120);norm(nMSEds_120);norm(nMSEl_120)];
%%
disp(biases_60');
disp(biases_90');
disp(biases_120');
disp(covariances_60');
disp(covariances_90');
disp(covariances_120');

Omega = [60 90 120];
biases = [biases_60 biases_90 biases_120];
covars = [covariances_60 covariances_90 covariances_120];
mses = [mse_60 mse_90 mse_120];
fig=figure('Renderer', 'painters', 'Position', [10 10 900 800])
f = subplot(2,2,1);
hold on
for i=1:size(biases,1)
    plot(Omega,biases(i,:));
end
title("Norm-Bias versus System Size ($\|\theta - \sum_{i=1}^N\hat{\theta_i}/N\|$ versus $\Omega$)",'Interpreter','latex','Fontsize',14 )
xlabel('$\Omega$','Interpreter','latex');
ylabel('$\|\theta - \sum_{i=1}^N\hat{\theta_i}/N\|$','Interpreter','latex');
hold off
g = subplot(2,2,2);
hold on
for i=1:size(covars,1)
    plot(Omega,log10(covars(i,:)));
end
title("Log-Covariance versus System Size ($\log_{10}(|\textrm{Cov}(\hat\theta)|)$ versus $\Omega$)",'Interpreter','latex','Fontsize',14 )
xlabel('$\Omega$','Interpreter','latex');
ylabel('$\log_{10}(|\textrm{Cov}(\hat\theta)|)$','Interpreter','latex');
legend('Norm D','Norm Ds', 'Norm L', 'SSA D', 'SSA Ds', 'SSA L');
hold off
h=subplot(2,2,[3,4]);
hold on
for i=1:size(mses,1)
    plot(Omega,log10(mses(i,:)));
end
title("Mean-Squared Error versus System Size ($\sum_{i=1}^N\|\theta-\hat\theta_i\|^2/N$ versus $\Omega$)",'Interpreter','latex','Fontsize',14 )
xlabel('$\Omega$','Interpreter','latex');
ylabel('$\sum_{i=1}^N\|\theta-\hat\theta_i\|^2/N$','Interpreter','latex');
hold off
saveas(fig,'500x8 Comparison Plots.png');
