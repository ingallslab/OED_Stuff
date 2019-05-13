addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all


%%variables
x_sym = SX.sym('x');
u_sym = SX.sym('u');
a_sym = SX.sym('a');
K_sym = SX.sym('K');
n_sym = SX.sym('n');
theta = [a_sym,K_sym,n_sym];

g = a_sym*(u_sym.^n_sym)./(K_sym.^n_sym+u_sym.^n_sym)-x_sym;

xstar = a_sym *(u_sym.^n_sym)./(u_sym.^n_sym+K_sym.^n_sym); 
%casadi variable for steady state

xStar = Function('xStar',{u_sym,theta},{xstar}); 
%tell matlab to consider xStar as a function with value xstar

uVals = linspace(0,40,100);
sVals = xStar(uVals,[a,K,n]); %s stands for steady state


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

objective_1 = @(th) -computeLikelihood(expData,inputs,th,deviations);
objective_2 = @(th) -compLikelihood_nocasadi(expData,inputs,th,deviations);
%options = optimset('PlotFcns',@optimplotfval);
options = optimset('Display','final');
min_1 = fminsearch(objective_1,[3.5,9.5,3.5],options);
min_2 = fminsearch(objective_2,[3.5,9.5,3.5],options);
disp(min_1);
disp(min_2);
disp(min_1 - theta);

figure;hold on
%   error bars
h1=fill([inputs';flipud(inputs')],[sparse(xStar(inputs,[a,K,n]))'-deviations';
    flipud(sparse(xStar(inputs,[a,K,n]))'+deviations')],[.4 .4 .7],'linestyle','none');%,'-','LineWidth',2);
set(h1,'facealpha',.3)
plot(uVals,sparse(sVals)); 
plot(uVals,sparse(xStar(uVals,min_1)),'b');
plot(uVals,sparse(xStar(uVals,min_2)),'k');
hold off

