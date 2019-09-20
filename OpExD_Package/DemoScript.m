
%these clear matlab memory and open windows 
close all
clear all

%% Adds Casadi to your path
addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.4')
import casadi.*

%% Set up the model in our experimental design library
%declare the number of parameters
theta=SX.sym('theta',4); 
%declare number of inputs
u=SX.sym('u',1);

%Define a simple Hill function, plus leak, as the model for the mean
%response
mean = theta(1) + (theta(2)*u(1).^theta(4)) ./(theta(3).^theta(4)+u(1).^theta(4));
%Define variance as constant (doesn't matter what it is because its
%constant)
var = 1000^2;

%Together the mean, var and distribution type define an observation
%variable
ObservedVar=struct('Type',{'Normal'}, 'Mean',{mean}, 'Covariance',{var});
%group all observation variables (sometimes more than one) together in a
%list
ObsVarList={ObservedVar};

%instantiate class
myModel=Model(ObsVarList,theta,u);

%% Fit the model

green_perc1=[1e-3 1 5 10 100];
mean_GFP1=[1429.3585, 1752.9652, 6101.2478, 12562.0651, 22504.4872];
std_GFP1=[5257.9221 5199.3618 6267.9641 9830.5494 14269.9423];

green_perc2=[1e-3 1 3 4 5 20 40 100];
mean_GFP2=[2651.2863 1807.8593 3345.9691 3937.829 4835.2193 8698.3805 13277.7421 14638.9662];
std_GFP2=[13993.5107	4882.2014	4531.9329	4082.4989	6304.1961	8404.4679 12915.8615 13496.0732];

yvals=[mean_GFP1 mean_GFP2];
yList={yvals};
uvals=[green_perc1 green_perc2];

theta_guess=[2000,10000,5,2];
lbd=[0 0 0 0];
ubd=[5000 20000 100 3];
theta_val=myModel.MaxLikFit(theta_guess,uvals,yList,lbd,ubd);

u_vals=[0:1:100];
y_vals=myModel.evalMu(theta_val, u_vals);

[X, Y] = meshgrid(logspace(0,2,200), linspace(0,1.2*theta_val(2),200));
Z=zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j)=myModel.evalPDF(Y(i,j),theta_val, X(i,j));
    end
end


%% Design Expeirments

fim0=zeros(length(theta_val));
for i=1:length(uvals)
    fim0=fim0+myModel.evalFIM(theta_val, uvals(i))
end

u_min=1;
u_max=100;
u_grid=[u_min:1:u_max];
NumPts=10;

D_approx=myModel.DesignApprox(theta_val, u_grid,'D',fim0);
D_exact=myModel.DesignExact(theta_val,NumPts,linspace(u_min,u_max,NumPts),[u_min; u_max],'D',fim0);

D_approx_noFim=myModel.DesignApprox(theta_val, u_grid,'D');
D_exact_noFim=myModel.DesignExact(theta_val,NumPts,linspace(u_min,u_max,NumPts),[u_min; u_max],'D');


%% Plot to summarize results

figure
subplot(7,1,1:3)
hold on
surf(X,Y,Z-max(max(Z)))
view(2)
shading interp
plot(u_vals,y_vals,'--k','LineWidth',3)
plot(uvals,yvals,'r*')
hold off
%myModel.evalLogLik(y_val,theta_val, u_val)

subplot(8,1,4)
stem(D_approx.uvals,D_approx.weights)

[us,cnts]=unique(round(D_exact.uvals,1));
subplot(8,1,5)
stem(us,cnts)

subplot(8,1,6)
stem(D_approx_noFim.uvals,D_approx_noFim.weights)

[us,cnts]=unique(round(D_exact_noFim.uvals,1));
subplot(8,1,7)
stem(us,cnts)

%add profile and asymptotic interval plots
%make everything easier to iterate
%add some plottin functions??

