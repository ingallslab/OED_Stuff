addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT COMPLETE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ll=logLik([2800.1261  86.7765],.5,theta,[.5 150],[phi_1;phi_2])

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

gamma=[.5 150];

pars=[theta gamma];

Omega=.01;

phi_1=[0.45 0];
phi_2=[0 45];

Phi=[phi_1;phi_2];

%% Casadi Stuff
x1_sym = SX.sym('x1_sym');
x2_sym = SX.sym('x2_sym');
x_sym=[x1_sym x2_sym];

% u1_sym = SX.sym('u1_sym');
% u2_sym = SX.sym('u2_sym');
% u_sym=[u1_sym u2_sym];
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

theta1_sym=[alpha1_sym beta1_sym K1_sym n1_sym kappa1_sym m1_sym];
theta2_sym=[alpha2_sym beta2_sym K2_sym n2_sym kappa2_sym m2_sym];
theta_sym=[alpha1_sym alpha2_sym beta1_sym beta2_sym K1_sym K2_sym n1_sym n2_sym kappa1_sym kappa2_sym m1_sym m2_sym];

F1 = alpha1_sym + beta1_sym./(1+(((alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym))./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;
g=[g1 g2];

F1_theta=jacobian(F1,theta_sym);
F1_x1=jacobian(F1,x1_sym);
x1_theta=-F1_theta./F1_x1;

g2_x1=jacobian(g2,x1_sym);
g2_theta=jacobian(g2,theta_sym);
x2_theta=g2_theta+g2_x1*x1_theta;

x_theta=[x1_theta;x2_theta];
%x_th_func = Function('x_th_func', {x_sym,u1_sym,u2_sym,theta_sym}, {x_theta});

b1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)+x1_sym;
b2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)+x2_sym;
B = [b1 0; 0 b2];
A=jacobian(g,x_sym);
A_func = Function('A_func',{x_sym,u_cntrl_sym,theta_sym},{A});

C_vec=inv(kron(eye(size(A)),A)+kron(A,eye(size(A))))*vec(-B);
C=reshape(C_vec,size(A));
C_func = Function('C_func', {x_sym,u_cntrl_sym,theta_sym}, {C});

%FIM=(x_theta.*repmat(theta_sym,2,1))'*C*(x_theta.*repmat(theta_sym,2,1));
FIM=(x_theta)'*C*(x_theta);
FIM_func = Function('FIM_func', {x_sym,u_sym,theta_sym}, {FIM});

x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
[x1_low,~,x1_high]=fixed_point_v5(x1_root,3000);
x_low=[x1_low x2_eval(x1_low)];
x_high=[x1_high x2_eval(x1_high)];

x1_low_sym = SX.sym('x1_low_sym');
x2_low_sym = SX.sym('x2_low_sym');
x_low_sym=[x1_low_sym x2_low_sym];

x1_hgh_sym = SX.sym('x1_hgh_sym');
x2_hgh_sym = SX.sym('x2_hgh_sym');
x_hgh_sym=[x1_hgh_sym x2_hgh_sym];

x1_obs_sym = SX.sym('x1_obs_sym');
x2_obs_sym = SX.sym('x2_obs_sym');
x_obs_sym=[x1_obs_sym x2_obs_sym];

gamma0_sym = SX.sym('gamma0_sym');
gamma1_sym = SX.sym('gamma1_sym');
gamma_sym=[gamma0_sym gamma1_sym];

%C_low_sym=C_func(x_low_sym,u_cntrl_sym,theta_sym);
%C_high_sym=C_func(x_hgh_sym,u_cntrl_sym,theta_sym);

%pdf_low=exp((x_obs_sym-x_low_sym)*inv(C_low_sym)*(x_obs_sym-x_low_sym)')/sqrt(det(C_low_sym)*(2*pi)^2);
%pdf_high=exp((x_obs_sym-x_hgh_sym)*inv(C_high_sym)*(x_obs_sym-x_hgh_sym)')/sqrt(det(C_high_sym)*(2*pi)^2);

pdf=exp(-.5*(x_obs_sym-x_sym)*inv(C)*(x_obs_sym-x_sym)')/sqrt(det(C)*(2*pi)^2);
pdf_func = Function('pdf_func', {x_obs_sym,x_sym,u_cntrl_sym,theta_sym}, {pdf});

pdf_low=pdf_func(x_obs_sym,x_low_sym,u_cntrl_sym,theta_sym);
pdf_high=pdf_func(x_obs_sym,x_hgh_sym,u_cntrl_sym,theta_sym);

u=0.45;
%u_val=(phi_2-phi_1)*u+phi_1;
C_low=full(C_func(x_low,u,theta))
C_high=full(C_func(x_high,u,theta))

mvnpdf([85.1261  986.7765],x_low,C_low)
pdf_func([85.1261  986.7765],x_low,u,theta)
mvnpdf([85.1261  986.7765],x_high,C_high)
pdf_func([85.1261  986.7765],x_high,u,theta)

rho=1/(1+exp(-(gamma1_sym*(u_cntrl_sym-gamma0_sym))));
ll=log(rho*pdf_low+(1-rho)*pdf_high);
ll_func = Function('ll_func', {x_obs_sym,x_low_sym,x_hgh_sym,u_cntrl_sym,theta_sym,gamma_sym}, {ll});

%logLik([85.1261  986.7765],u,theta,[.5 150],[phi_1;phi_2])
ll_func([85.1261  986.7765],x_low,x_high,u,theta,[.5 150])

% I_low=full(FIM_func(x_low,u_vals,theta));
% I_hgh=full(FIM_func(x_high,u_vals,theta));
% 
% cond(I_low)
% cond(I_hgh)
u_vec=[.1 .45 .9];
bounds=[.4 .5];
evalFIM(pars,Omega,u_vec,Phi,bounds)


function FIM_vec=evalFIM(pars,Omega,u_vec,Phi,bounds)

    addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*

    phi_1=Phi(1,:);
    phi_2=Phi(2,:);
    
    theta=pars(1:12);
    
    alpha_1 = theta(1);
    alpha_2 = theta(2); 
    beta_1 = theta(3);   
    beta_2 = theta(4);   
    K_1 = theta(5);   
    K_2 = theta(6);    
    n_1  = theta(7);    
    n_2 = theta(8);    
    kappa_1 = theta(9); 
    kappa_2 = theta(10);  
    m_1= theta(11);    
    m_2 = theta(12);
    
    gamma=pars(13:end);
       
    N_u=length(u_vec);
    Np=14;
    N_c=2;
    N_th=12;

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
    
    F1 = alpha1_sym + beta1_sym./(1+(((alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym))./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
    g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
    g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;
    g=[g1 g2];

    F1_theta=jacobian(F1,theta_sym);
    F1_x1=jacobian(F1,x1_sym);
    x1_theta=-F1_theta./F1_x1;

    g2_x1=jacobian(g2,x1_sym);
    g2_theta=jacobian(g2,theta_sym);
    x2_theta=g2_theta+g2_x1*x1_theta;

    x_theta=[x1_theta;x2_theta];
    x_theta_func = Function('x_theta_func', {x_sym,u_cntrl_sym,theta_sym}, {x_theta});

    b1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)+x1_sym;
    b2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)+x2_sym;
    B = [b1 0; 0 b2];
    A=jacobian(g,x_sym);
    %?????????????????????????????????????????????????????????
    C_vec=(inv(kron(eye(size(A)),A)+kron(A,eye(size(A))))*vec(-B))./Omega_sym;
    %?????????????????????????????????????????????????????????
    C=reshape(C_vec,size(A)); 
    Cinv=inv(C);
    Cinv_func = Function('Cinv_func', {x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {Cinv});
    C_theta=jacobian(C_vec([1 2 4]),theta_sym);
    %FIM=(x_theta.*repmat(theta_sym,2,1))'*C*(x_theta.*repmat(theta_sym,2,1));
    
    F=[x_theta;C_theta];
    D=[1 0 0; 0 1 0; 0 1 0; 0 0 1];
    FIM=F'* [inv(C) zeros(2,3); zeros(3,2)  .5*D'*kron(inv(C),inv(C))*D]*F ;%(x_theta)'*inv(C)*(x_theta);
    FIM_func = Function('FIM_func', {x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {FIM});


    pdf=exp(-.5*(x_obs_sym-x_sym)*inv(C)*(x_obs_sym-x_sym)')/sqrt(det(C)*(2*pi)^2);
    pdf_func = Function('pdf_func', {x_obs_sym,x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {pdf});

    pdf_low=pdf_func(x_obs_sym,x_low_sym,u_cntrl_sym,theta_sym,Omega_sym);
    pdf_high=pdf_func(x_obs_sym,x_high_sym,u_cntrl_sym,theta_sym,Omega_sym);


    rho=1/(1+exp(-(gamma1_sym*(u_cntrl_sym-gamma0_sym))));
    Lik=rho*pdf_low+(1-rho)*pdf_high;
    Lik_func = Function('Lik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {Lik});
    
    logLik=log(Lik);
    logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {logLik});

    %log sensitivities
    logLik_pars=jacobian(logLik,par_sym);
    %logLik_pars=jacobian(logLik,par_sym).*par_sym;

    logLik_xlow=jacobian(logLik,x_low_sym);
    logLik_xhigh=jacobian(logLik,x_high_sym);

    logLik_sens=logLik_pars+[logLik_xlow*x_theta_func(x_low_sym,u_cntrl_sym,theta_sym) zeros(1,N_c)]+[logLik_xhigh*x_theta_func(x_high_sym,u_cntrl_sym,theta_sym) zeros(1,N_c)];

    logLik_sens_func = Function('logLik_sens_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {logLik_sens});

    FIM_vec=zeros(Np,Np,N_u);

    for i=1:N_u

        if u_vec(i)<=bounds(1)
            %in lower monostable region
            u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
            x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
            x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
            [x1_low,~,~]=fixed_point_v5(x1_root,3000);
            x_low=[x1_low x2_eval(x1_low)];
            %x_high=[x1_high x2_eval(x1_high)];
            

            %FIM_vec(:,:,i)=(x_th'*x_th)./sig2+0.5*trace(sig2_th'*sig2_th)./(sig2^2);
            FIM_vec(1:N_th,1:N_th,i)=full(FIM_func(x_low,u_vec(i),theta,Omega));
            %(x_th'*x_th)./sig2+0.5*((sig2_th+sig2_xval.*x_th)./sig2)'*((sig2_th+sig2_xval.*x_th)./sig2);

        elseif (bounds(1)<u_vec(i))&&(u_vec(i)<bounds(2))
            %in bistable region
            u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
            x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
            x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
            [x1_low,~,x1_high]=fixed_point_v5(x1_root,3000);
            x_low=[x1_low x2_eval(x1_low)];
            x_high=[x1_high x2_eval(x1_high)];
            
            %Cinv_low = Cinv_func(x_sym,u_cntrl_sym,theta_sym,Omega_sym);
            %Cinv_high = Cinv_func(x_sym,u_cntrl_sym,theta_sym,Omega_sym);
            
            %std_low=sqrt(full(sigma2_func(low,u_vec(i),theta,Omega)));
            %std_hgh=sqrt(full(sigma2_func(high,u_vec(i),theta,Omega)));
            
            lB1=0;%min(low-4*std_low,high-4*std_hgh);
            uB1=2000;%max(low+4*std_low,high+4*std_hgh);
            lB2=0;
            uB2=2000;

            FIM_integrand=@(x1,x2)...
                (full(logLik_sens_func([x1,x2],x_low,x_high,u_vec(i),pars,Omega))'...
                *full(logLik_sens_func([x1,x2],x_low,x_high,u_vec(i),pars,Omega)))...
                *full(Lik_func([x1,x2],x_low,x_high,u_vec(i),pars,Omega));
            %FIM_full = integral(FIM_integrand,-2e1,2e1,'ArrayValued',true);
            %should make bounds a function of the stnd dev,omega etc.
            FIM_full = integral2(FIM_integrand,lB1,uB1,lB2,uB2);

            FIM_vec(:,:,i)=FIM_full;

        else
            %in upper monostable region
            u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
            x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
            x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
            [~,~,x1_high]=fixed_point_v5(x1_root,3000);
            x_high=[x1_high x2_eval(x1_high)];
            %x_high=[x1_high x2_eval(x1_high)];
            

            %FIM_vec(:,:,i)=(x_th'*x_th)./sig2+0.5*trace(sig2_th'*sig2_th)./(sig2^2);
            FIM_vec(1:N_th,1:N_th,i)=full(FIM_func(x_high,u_vec(i),theta,Omega));
            %(x_th'*x_th)./sig2+0.5*((sig2_th+sig2_xval.*x_th)./sig2)'*((sig2_th+sig2_xval.*x_th)./sig2);

        end
        FIM_vec(:,:,i)=FIM_vec(:,:,i).*(pars'*pars);
    end
end

% function ll=logLik(x,u,theta,gamma,Phi)
% 
%     addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
%     import casadi.*
% 
%     x1_sym = SX.sym('x1_sym');
%     x2_sym = SX.sym('x2_sym');
%     x_sym=[x1_sym x2_sym];
% 
%     u1_sym = SX.sym('u1_sym');
%     u2_sym = SX.sym('u2_sym');
%     
%     u_sym=[u1_sym u2_sym];
%     Omega_sym = SX.sym('Omega_sym');
% 
%     alpha1_sym = SX.sym('alpha1_sym');
%     alpha2_sym = SX.sym('alpha2_sym'); 
%     beta1_sym = SX.sym('beta1_sym');  
%     beta2_sym = SX.sym('beta2_sym');     
%     K1_sym = SX.sym('K1_sym');    
%     K2_sym = SX.sym('K2_sym');   
%     n1_sym  = SX.sym('n1_sym');     
%     n2_sym = SX.sym('n2_sym');    
%     kappa1_sym = SX.sym('kappa1_sym');  
%     kappa2_sym = SX.sym('kappa2_sym');   
%     m1_sym = SX.sym('m1_sym');     
%     m2_sym = SX.sym('m2_sym');
% 
%     theta1_sym=[alpha1_sym beta1_sym K1_sym n1_sym kappa1_sym m1_sym];
%     theta2_sym=[alpha2_sym beta2_sym K2_sym n2_sym kappa2_sym m2_sym];
%     theta_sym=[alpha1_sym alpha2_sym beta1_sym beta2_sym K1_sym K2_sym n1_sym n2_sym kappa1_sym kappa2_sym m1_sym m2_sym];
% 
%     g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u2_sym./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
%     g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u1_sym./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;
% 
%     g=[g1 g2];
%     g_x=jacobian(g,x_sym);
%     g_x_func = Function('g_x_func', {x_sym,u_sym,theta_sym}, {g_x});
% 
%     alpha_1 = theta(1);
%     alpha_2 =theta(2);
%     beta_1 = theta(3);
%     beta_2 = theta(4);
%     K_1 = theta(5);
%     K_2 = theta(6);
%     n_1  = theta(7);
%     n_2 = theta(8);
%     kappa_1 = theta(9); 
%     kappa_2 = theta(10); 
%     m_1= theta(11);     
%     m_2 = theta(12);
%     
%     phi_1=Phi(1,:);
%     phi_2=Phi(2,:);
%     u_vals=(phi_2-phi_1)*u+phi_1;
% 
%     gamma_0=gamma(1);
%     gamma_1=gamma(2);
%     
%     x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
%     x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
%     B_func=@(x) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)+x(1) 0;...
%             0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2)+x(2)];
%     
%     [x1_low,~,x1_high]=fixed_point_v5(x1_root,3000);
%     x_low=[x1_low x2_eval(x1_low)];
%     x_high=[x1_high x2_eval(x1_high)];
%     
%     A_low=full(g_x_func(x_low,u_vals,theta));
%     B_low=B_func(x_low);
%     C_low = lyap(A_low,B_low);
%     
%     A_hgh=full(g_x_func(x_high,u_vals,theta));
%     B_hgh=B_func(x_high);
%     C_hgh = lyap(A_hgh,B_hgh);
%     
%     rho=1/(1+exp(-(gamma_1*(u-gamma_0))));
% 
%     ll=log(rho*mvnpdf(x,x_low,C_low)+(1-rho)*mvnpdf(x,x_high,C_hgh));
% 
% end
