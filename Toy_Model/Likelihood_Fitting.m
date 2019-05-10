addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

expData = dlmread('./Toy_Model/Data/SSA_Data_90.txt','\t');
inputs = dlmread('./Toy_Model/Data/SSA_Data_90_uVals.txt','\t');

deviations=linspace(1,size(expData,2),size(expData,2));
for i=1:size(expData,2)
    deviations(i) = std(expData(:,i));
end

a=3;
K=9;
n=3;
theta = [a,K,n];

disp(computeLikelihood(expData,inputs,[a,K,n],deviations));
disp(compLikelihood_nocasadi(expData,inputs,theta,deviations));

objective = @(th) -compLikelihood_nocasadi(expData,inputs,th,deviations);
options = optimset('PlotFcns',@optimplotfval);
min = fminsearch(objective,[2.5,9.5,3.2],options);

disp(min);

disp(min - theta);