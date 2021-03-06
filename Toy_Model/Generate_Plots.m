addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%some nice values
a=3;
K=9;
n=3;
m=2;

%Casadi Setup
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

%%
%%%-----------------Old Plot Generator-----------------

%zVals = xStar(uVals,[a,K,m]); %I wanted to see what n=2 looks like
%
%%sparse converts casadi.DM classtype to a list of doubles
%
%generating normally distributed data points (~N(0,0.1))
%normRands = normrnd(0,0.1,[100,1]);
%Redist_sVals = abs(normRands+sVals(:));%normally redistributed sVals
%scatter(uVals,sparse(Redist_sVals),'filled');
%
%%generating more interesting normally distributed data points
%xNormRands = linspace(0,1,100); %preallocating array for efficiency
%%x dependent random variables -> N(0,x*0.2)
%for i=1:100
%   xNormRands(i) = normrnd(0,0.2*sVals(i));
%end
%xRedist_sVals = abs(xNormRands+sVals(:));
%scatter(uVals,sparse(xRedist_sVals),'filled');
%
%%properly generating normally distributed "experimental" data
%
%inputVals_p = linspace(0,40,15);
%outputVals_p = sparse(xStar(inputVals_p,[a,K,n]));
%
%properVals_p = meshgrid(1:15,1:15);
%for i=1:15
%    properVals_p(:,i) = abs(outputVals_p(i)+normrnd(0,outputVals_p(i)*0.08,[15,1]));
%end
%hold on
%for i=1:15
%    scatter(ones(1,15)*inputVals_p(i),properVals_p(:,i));
%end
%hold off
%%
%%%-----------------Stochastic Simulation Algorithm-----------------

numExperiments = 40;
numSamples = 20;

generatedData = meshgrid(1:numExperiments,1:numSamples);
inputVals = linspace(0,40,numExperiments);

Omega = 90;
finTime = 1500;

for i=1:numExperiments
    generatedData(:,i)=SSA_Func(zeros(numSamples,1),inputVals(i),Omega,finTime,numSamples);
end

deviations=linspace(1,numExperiments,numExperiments);
for i=1:numExperiments
    deviations(i) = std(generatedData(:,i));
end

%disp(generatedData);
%for i=1:10
    %scatter(ones(40,1)*inputVals(i),generatedData(:,i));
%end
figure;hold on
%   error bars
h1=fill([inputVals';flipud(inputVals')],[sparse(xStar(inputVals,[a,K,n]))'-deviations';
    flipud(sparse(xStar(inputVals,[a,K,n]))'+deviations')],[.4 .4 .7],'linestyle','none');%,'-','LineWidth',2);
set(h1,'facealpha',.3)
plot(uVals,sparse(sVals)); 
hold off

dlmwrite('./Toy_Model/Data/SSA_Data_90.txt',generatedData,'\t');
dlmwrite('./Toy_Model/Data/SSA_Data_90_uVals.txt',inputVals,'\t');










