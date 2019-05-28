numExperiments = 128;
numSamples = 100;

generatedData = meshgrid(1:numExperiments,1:numSamples);
inputVals = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];

Omega = 300;
finTime = 12000;

for i=1:16:numExperiments
    generatedData(:,i)=SSA_Generator_1D(zeros(numSamples,1),inputVals(i),Omega,finTime,numSamples);
    disp(i);
end

% dlmwrite('./1D_Michael/Data/SSA_Data_300.txt',generatedData,'\t');
% dlmwrite('./1D_Michael/Data/SSA_Data_300_uVals.txt',inputVals,'\t');

hold on
plotbif(u_vec,theta_t);
for i=1:numExperiments
    scatter(inputVals(i)*ones(numSamples,1),generatedData(:,i));
end
hold off