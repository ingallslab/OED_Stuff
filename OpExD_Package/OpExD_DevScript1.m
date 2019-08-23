close all
clear all
addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*



%% MODEL 1
% MODEL 1 DEF: 
% [y1,y2]~Norm([mu11,mu12],[mu13, mu15; mu15, mu14])
% y3~Bern(mu2)
% y4~Norm(mu31,mu32)

%set inputs and parameter vecs
theta=SX.sym('theta',7); %define model parameters
u=SX.sym('u',2); %define experimental inputs

%define model for pdf params

% Example for two models: 
% models={ struct('Type',{'Normal'}, 'Mean',{mu1_sym}, ...
%           'Covariance',{cov1_sym}, 'Details',{'2x2 normpdf'},...
%           'Params',{theta1_sym}, 'Input',{u1_sym}),
%          struct('Type',{'LogNormal'}, 'Mean',{mu2_sym}, ...
%           'Covariance',{cov2_sym}, 'Details',{'4z4 lognormpdf'},...
%           'Params',{theta2_sym}, 'Input',{u2_sym})}
%Here, sym_vec is a vector containing the individual symbols
%required for the models - this allows us to build the loglik_func 
%and fim_func objects. Details is a string describing the model for debugging purposes.
%Params are the desired parameters we need to fit to.
                
mu11 = theta(1) + theta(2)*u(1) +  theta(3)*u(2) + theta(4)*u(1)*u(2);
mu12 = theta(5) + theta(6)*u(1) +  theta(7)*u(2);
mu13 = theta(2)*u(1)^2;
mu14 = theta(3)*u(2)^2;
mu15 = theta(4)*u(1)*u(2);
mu2 =  theta(1)/(1+theta(1));
mu31 = theta(1) + theta(6)*u(1) +  theta(7)*u(2);
mu32 = theta(1);
mean1=[mu11; mu12];
covar1 = [[mu13 mu15];[mu15 mu14]];
params1=struct('Type',{'Normal'}, 'Mean',{mean1}, 'Covariance',{covar1},...
               'Params',{theta}, 'Input',{u}, 'Details',{'normpdf 2x2'});
mean2=mu2;
params2=struct('Type',{'Normal'}, 'Mean',{mean2}, 'Covariance',{1},...
               'Params',{theta}, 'Input',{u}, 'Details',{'normpdf 1x1, default covariance'});
mean3 = theta(1) + theta(6)*u(1) +  theta(7)*u(2);
covar3 = theta(1);
params3=struct('Type',{'Normal'}, 'Mean',{mean3}, 'Covariance',{covar3},...
               'Params',{theta}, 'Input',{u}, 'Details',{'normpdf 1x1'});

%assemble model
model={params1, params2, params3};

%instantiate class
Model1=OpExD_Model(model,theta,u);

%MODE 1 TEST:
theta_val=[1,1,1,1,1,1,1]; u_val=[1,1]; x_val=[1,1];%nominal params - don't work atm
Model1.evalPDF(x_val,theta_val, u_val)
Model1.evalLogLik(x_val,theta_val, u_val)
Model1.evalFIM(theta_val, u_val)
u_grid=[0:0.1:3];
expD=Model1.optRelax(theta_val, u_grid,'D');
expD=Model1.optRelax(theta_val, u_grid,'Ds',2);
expA=Model1.optRelax(theta_val, u_grid,'A');
expE=Model1.optRelax(theta_val, u_grid,'E');
expV=Model1.optRelax(theta_val, u_grid,'V',[0.1:0.1:0.3]);
expV=Model1.optRelax(theta_val, u_grid,'V',[3:0.1:3.5]);

expD_exact=expD.genExact(10);

exactD=Model.optExact(theta_val,10,linspace(0,3,10),[0; 3],'D')
exactA=Model.optExact(theta_val,10,linspace(0,3,10),[0; 3],'A')


clear theta u mu Model theta_val u_val x_val
% %% MODEL 2
% 
% % MODEL 2: Hill function, single input, proportional var (with known const)
% theta=SX.sym('theta',3); %define model parameters
% u=SX.sym('u'); %define experimental inputs
% mu=[theta(1)*u/(theta(2)^theta(3)+u^theta(3))...
%     .1*(theta(1)*u/(theta(2)^theta(3)+u^theta(3)))]; %define model for mean and varaince 
% Model=OpExD_Model(theta,u,mu,'Norm1_MeanVar');
% 
% %MODE 1 TEST
% theta_val=[3 .5 2]; u_val=2; x_val=1;%nominal params
% Model.evalPDF(x_val,theta_val, u_val)
% Model.evalLogLik(x_val,theta_val, u_val)
% Model.evalFIM(theta_val, u_val)
% u_grid=[0.1:0.1:3];
% expD=Model.optRelax(theta_val, u_grid,'D');
% expA=Model.optRelax(theta_val, u_grid,'A');
% expE=Model.optRelax(theta_val, u_grid,'E');
% expV=Model.optRelax(theta_val, u_grid,'V',[0.1:0.1:0.3]);
% 
% exactD=Model.optExact(theta_val,3,linspace(0.1,3,3),[0.01; 5],'D')
% exactA=Model.optExact(theta_val,3,linspace(0.1,3,3),[0.01; 5],'A')
% 
% clear theta u mu Model theta_val u_val x_val
%% MODEL 3




