addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

expData = dlmread('./Toy_Model/Data/SSA_Data_90.txt','\t');

%casadi setup
x_sym = SX.sym('x');
u_sym = SX.sym('u');
a_sym = SX.sym('a');
K_sym = SX.sym('K');
n_sym = SX.sym('n');
theta = [a_sym,K_sym,n_sym];
sigma_sym = SX.sym('sigma');

xstar = a_sym *(u_sym.^n_sym)./(u_sym.^n_sym+K_sym.^n_sym); 

x_obs_sym = SX.sym('x_obs');

lik = (1/(sqrt(2*pi)*sigma_sym))*exp(-((x_obs_sym-xstar).^2)./(2*sigma_sym.^2));
loglik = log(lik);
loglik_f = Function('loglik_f',{x_obs_sym,u_sym,theta_sym,sigma_sym},{loglik});

loglik_th = jacobian(loglik,theta_sym);
loglik_th_f = Function('loglik_theta_f',{x_obs_sym,u_sym,theta_sym,sigma_sym},{loglik_th});

deviations=linspace(1,30,30);
for i=1:30
    deviations(i) = std(expData(:,i));
end

