close all
clear all
addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*



%% MODEL 1
% MODEL 1 DEF: 
% [y1,y2]~Norm([mu11,mu12],[mu13, mu15; mu15, mu14])
% y3~Bern(mu2)
% y4~Norm(mu31,mu32)
%mu11 = theta(1) + theta(2)*u(1) +  theta(3)*u(2) + theta(4)*u(1)*u(2)
%mu12 = theta(5) + theta(6)*u(1) +  theta(7)*u(2)
%mu13 = theta(2)*u(1)^2
%mu14 = theta(3)*u(2)^2
%mu15 = theta(4)*u(1)*u(2)
%mu2 =  theta(1)/(1+theta(1))
%mu31 = theta(1) + theta(6)*u(1) +  theta(7)*u(2)
%mu32 = theta(1) 

%set inputs and parameter vecs
theta=SX.sym('theta',7); %define model parameters
u=SX.sym('u',2); %define experimental inputs

%define model for pdf params
mu1=[theta(1) + theta(2)*u(1) +  theta(3)*u(2) + theta(4)*u(1)*u(2),... $mean 1
    theta(5) + theta(6)*u(1) +  theta(7)*u(2),... %mean 2
    theta(2)*u(1)^2, theta(3)*u(2)^2,theta(4)*u(1)*u(2)]; %covariance parts
mu2=[theta(1)/(1+theta(1))]; 
mu3=[theta(1) + theta(6)*u(1) +  theta(7)*u(2), theta(1) ];

%assemble model
model={{'Norm',mu1},{'Norm',mu2},{'Norm',mu3}};
%instantiat class
Model1=OpExD_Model(model,theta,u);



%MODE 1 TEST:
theta_val=[4 2 .5]; u_val=2; x_val=1;%nominal params
Model.evalPDF(x_val,theta_val, u_val)
Model.evalLogLik(x_val,theta_val, u_val)
Model.evalFIM(theta_val, u_val)
u_grid=[0:0.1:3];
expD=Model.optRelax(theta_val, u_grid,'D');
expD=Model.optRelax(theta_val, u_grid,'Ds',2);
expA=Model.optRelax(theta_val, u_grid,'A');
expE=Model.optRelax(theta_val, u_grid,'E');
expV=Model.optRelax(theta_val, u_grid,'V',[0.1:0.1:0.3]);
expV=Model.optRelax(theta_val, u_grid,'V',[3:0.1:3.5]);

expD_exact=expD.genExact(10);

exactD=Model.optExact(theta_val,10,linspace(0,3,10),[0; 3],'D')
exactA=Model.optExact(theta_val,10,linspace(0,3,10),[0; 3],'A')


clear theta u mu Model theta_val u_val x_val
%% MODEL 2

% MODEL 2: Hill function, single input, proportional var (with known const)
theta=SX.sym('theta',3); %define model parameters
u=SX.sym('u'); %define experimental inputs
mu=[theta(1)*u/(theta(2)^theta(3)+u^theta(3))...
    .1*(theta(1)*u/(theta(2)^theta(3)+u^theta(3)))]; %define model for mean and varaince 
Model=OpExD_Model(theta,u,mu,'Norm1_MeanVar');

%MODE 1 TEST
theta_val=[3 .5 2]; u_val=2; x_val=1;%nominal params
Model.evalPDF(x_val,theta_val, u_val)
Model.evalLogLik(x_val,theta_val, u_val)
Model.evalFIM(theta_val, u_val)
u_grid=[0.1:0.1:3];
expD=Model.optRelax(theta_val, u_grid,'D');
expA=Model.optRelax(theta_val, u_grid,'A');
expE=Model.optRelax(theta_val, u_grid,'E');
expV=Model.optRelax(theta_val, u_grid,'V',[0.1:0.1:0.3]);

exactD=Model.optExact(theta_val,3,linspace(0.1,3,3),[0.01; 5],'D')
exactA=Model.optExact(theta_val,3,linspace(0.1,3,3),[0.01; 5],'A')

clear theta u mu Model theta_val u_val x_val
%% MODEL 3




