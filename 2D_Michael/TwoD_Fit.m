addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all  

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

u0_1=0.25;
u0_2=20;
x1_null=@(x_2,u_2) alpha_1 + beta_1./(1+((x_2./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1);
x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);

g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
g_func2=@(x_2,u_1,u_2) alpha_2 + beta_2*1./(1+(((alpha_1 +beta_1*1./(1+((x_2./K_2)*(1/(1+(u_2./kappa_2).^m_2))).^n_1))./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x_2;

g_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x(1) alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x(2)];

B_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)+x(1) 0;...
    0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)+x(2)];

umax_1=1;
umax_2=100;
[u1_grid,u2_grid] = meshgrid(0:0.1:umax_1,0:10:umax_2);
numTrials=20;
sy = TwoD_Symbols_EM(u1_grid,u2_grid,numTrials);
solver = nlpsol('solver','ipopt',sy.nlp);


SSAData=cell(size(u1_grid,1),size(u1_grid,2));

x1vals=[];
x2vals=[];
u1vals=[];
u2vals=[];

x0=theta';

for i=1:size(u1_grid,1)
    for j=1:size(u2_grid,2)
        s = dlmread(strcat('2D_Michael/Data/SSAData_',num2str(i),'_',num2str(j)),'\t');
        SSAData{i,j}=s;
        
        g_roots1=@(x_1) g_func1(x_1,u1_grid(i,j),u2_grid(i,j));
        [x1_low,~,x1_high]=fixed_point_v5(g_roots1,3000);
        x2_low=x2_null(x1_low,u1_grid(i,j));
        x2_high=x2_null(x1_high,u1_grid(i,j));
        
        x1vals=[x1vals; s(:,1)];
        x2vals=[x2vals; s(:,2)];
        u1vals=[u1vals;u1_grid(i,j)];
        u2vals=[u2vals;u2_grid(i,j)];
        
        x0=[x0; x1_high];
        x0=[x0; x2_high];
        x0=[x0; x1_low];
        x0=[x0; x2_low];
        
        x0=[x0; 1];  %  covariance matrix
        x0=[x0; 0];
        x0=[x0; 0];
        x0=[x0; 1];
        
        x0=[x0; 1];  %  covariance matrix
        x0=[x0; 0];
        x0=[x0; 0];
        x0=[x0; 1];
        
        x0=[x0;0.5]; %  distribution weighting
        x0=[x0;0.5];
    end
end
dataSet=[x1vals;x2vals;u1vals;u2vals];

%Fit=solver('x0',x0,'lbx',sy.lbx,'ubx',sy.ubx,'lbg',sy.lbg,'ubg',sy.ubg,'p',dataSet);

disp('test');