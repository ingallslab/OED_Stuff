numExperiments = 128;
numSamples = 40000;

generatedData = meshgrid(1:128,1:numSamples);

Omega = 120;
finTime = 8000;
u_vec = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
%87
parfor i=80:128
    generatedData(:,i)=SSA_Compute_1D(1.5*ones(numSamples,1),u_vec(i),Omega,finTime,numSamples);
    disp(i);
    dlmwrite(strcat('./1D_Michael/drive_W120/SSA_Data_120_u=',num2str(u_vec(i))),generatedData(:,i),'\t');
end

generatedData2 = meshgrid(1:numExperiments,1:numSamples);
Omega2 = 60;
finTime2 = 5000;
parfor i=1:128
    generatedData2(:,i)=SSA_Compute_1D(1.5*ones(numSamples,1),u_vec(i),Omega2,finTime2,numSamples);
    disp(i);
    dlmwrite(strcat('./1D_Michael/drive_W60/SSA_Data_60_u=',num2str(u_vec(i))),generatedData2(:,i),'\t');
end
