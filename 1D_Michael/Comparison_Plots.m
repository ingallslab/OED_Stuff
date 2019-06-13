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

mse_90 = [sMSEd_90;sMSEds_90;sMSEl_90;nMSEd_90;nMSEds_90;nMSEl_90];
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

mse_60 = [sMSEd_60;sMSEds_60;sMSEl_60;nMSEd_60;nMSEds_60;nMSEl_60];
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

mse_120 = [sMSEd_120;sMSEds_120;sMSEl_120;nMSEd_120;nMSEds_120;nMSEl_120];
%%
disp(biases_60');
disp(biases_90');
disp(biases_120');
disp(covariances_60');
disp(covariances_90');
disp(covariances_120');





