addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*

phi_1=[0.45 0];
phi_2=[0 45];
Phi=[phi_1;phi_2];
x1_sym = SX.sym('x1_sym');
x2_sym = SX.sym('x2_sym');
x_sym=[x1_sym x2_sym];
Omega_sym = SX.sym('Omega_sym');
u_cntrl_sym = SX.sym('u_cntrl_sym');
u_val_sym=(phi_2-phi_1)*u_cntrl_sym+phi_1;

alpha1_sym = SX.sym('alpha1_sym');
alpha2_sym = SX.sym('alpha2_sym'); 
beta1_sym = SX.sym('beta1_sym');  
beta2_sym = SX.sym('beta2_sym');     
K1_sym = SX.sym('K1_sym');    
K2_sym = SX.sym('K2_sym');   
n1_sym  = SX.sym('n1_sym');     
n2_sym = SX.sym('n2_sym');    
kappa1_sym = SX.sym('kappa1_sym');  
kappa2_sym = SX.sym('kappa2_sym');   
m1_sym = SX.sym('m1_sym');     
m2_sym = SX.sym('m2_sym');
theta_sym=[alpha1_sym alpha2_sym beta1_sym beta2_sym K1_sym K2_sym n1_sym n2_sym kappa1_sym kappa2_sym m1_sym m2_sym];

gamma0_sym = SX.sym('gamma0_sym');
gamma1_sym = SX.sym('gamma1_sym');
gamma_sym=[gamma0_sym gamma1_sym];

par_sym=[theta_sym gamma_sym];

x1_low_sym = SX.sym('x1_low_sym');
x2_low_sym = SX.sym('x2_low_sym');
x_low_sym=[x1_low_sym x2_low_sym];

x1_high_sym = SX.sym('x1_high_sym');
x2_high_sym = SX.sym('x2_high_sym');
x_high_sym=[x1_high_sym x2_high_sym];

x1_obs_sym = SX.sym('x1_obs_sym');
x2_obs_sym = SX.sym('x2_obs_sym');
x_obs_sym=[x1_obs_sym x2_obs_sym];

g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;
g=[g1 g2];

b1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)+x1_sym;
b2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)+x2_sym;
B = [b1 0; 0 b2];
A=jacobian(g,x_sym);
A_func = Function('A_func',{x_sym,u_cntrl_sym,theta_sym},{A});

B_func = Function('B_func',{x_sym,u_cntrl_sym,theta_sym},{B});

C_vec=(inv(kron(eye(size(A)),A)+kron(A,eye(size(A))))*vec(-B))./Omega_sym;
Cvec_func = Function('Cvec_func',{x_sym,u_cntrl_sym,theta_sym,Omega_sym},{C_vec});

C=reshape(C_vec,size(A)); 
Cinv=inv(C);
Cinv_func = Function('Cinv_func', {x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {Cinv});
C_theta=jacobian(C_vec([1 2 4]),theta_sym);

pdf=exp(-.5*(x_obs_sym-x_sym)*inv(C)*(x_obs_sym-x_sym)')/sqrt(det(C)*(2*pi)^2);
pdf_func = Function('pdf_func', {x_obs_sym,x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {pdf});

pdf_low=pdf_func(x_obs_sym,[x1_low_sym,x2_high_sym],u_cntrl_sym,theta_sym,Omega_sym);
pdf_high=pdf_func(x_obs_sym,[x1_high_sym,x2_low_sym],u_cntrl_sym,theta_sym,Omega_sym);

rho=1/(1+exp(-(gamma1_sym*(u_cntrl_sym-gamma0_sym))));
rho_func=Function('rho_func',{par_sym,u_cntrl_sym},{rho});
Lik=rho*pdf_low+(1-rho)*pdf_high;
Lik_func = Function('Lik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {Lik});

logLik=log(Lik);
logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {logLik});

logLik_xH = gradient(logLik,x_high_sym);
logLik_xH_func = Function('logLik_xH_func',{x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym},{logLik_xH});

uVals=0.05:0.1:0.95;
SSAData=cell(length(uVals),1);

syms=struct('pdf_func',pdf_func,'rho_func',rho_func);

for i=1:length(uVals)
    u1 = uVals(i);
    uV = (phi_2-phi_1)*u1+phi_1;
    u1 = uV(1);
    u2 = uV(2);
    s = dlmread(strcat('2D_Michael/Data/SliceData_Omega=0.01_u=',num2str(u1)),'\t');
    SSAData{i}=s(1:20,:);
end

Omega_true=0.01;
alpha_1 = 13.609;
alpha_2 = 60.882;   
beta_1 = 3529.923;   
beta_2 = 1053.916;      
K_1 = 31.94;    
K_2 = 30.0;    
n_1  = 2.00;     
n_2 = 2.00;     
kappa_1 = 0.0906;   
kappa_2 = 11.65;    
m_1= 2.00;     
m_2 = 2.00;

theta=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];
theta_true= [theta .5 150];

objec=@(theta) Objective(uVals,SSAData,theta,Omega_true,syms);
opti=optimset('Display','iter','FunValCheck','off','MaxFunEvals',814);
mini=fminsearch(objec,theta_true,opti);
disp(mini);

disp('pause');

function obje=Objective(uVals,xVals,theta,Omega,s)
    ll=0;
    phi_1=[0.45 0];
    phi_2=[0 45];
    Phi=[phi_1;phi_2];
    
    %theta=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];
    
    alpha_1=theta(1);
    alpha_2=theta(2);
    beta_1=theta(3);
    beta_2=theta(4);
    K_1=theta(5);
    K_2=theta(6);
    n_1=theta(7);
    n_2=theta(8);
    kappa_1=theta(9);
    kappa_2=theta(10);
    m_1=theta(11);
    m_2=theta(12);

    x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);
    g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
    
    for i=1:length(uVals)
        for j=1:size(xVals{i},2)
            xV=xVals{i}(j,:);
            u1 = uVals(i);
            uV = (phi_2-phi_1)*u1+phi_1;
            u1 = uV(1);
            u2 = uV(2);

            g_roots1=@(x_1) g_func1(x_1,u1,u2);
            [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);

            if (x1_low==x1_mid&&x1_mid==x1_high)    
                x2_low=x2_null(x1_low,u1);
                ll=ll+full(log(s.pdf_func(xV,[x1_low x2_low],u1,theta(1:12),Omega)));
            else
                x2_low=x2_null(x1_high,u1);
                x2_high=x2_null(x1_low,u1);

                ll=ll+full(log(s.rho_func(theta,u1)*s.pdf_func(xV,[x1_low x2_high],u1,theta(1:12),Omega)+...
                    (1-s.rho_func(theta,u1))*s.pdf_func(xV,[x1_high x2_low],u1,theta(1:12),Omega)));
            end
        end
    end
    obje=-ll;
end
    