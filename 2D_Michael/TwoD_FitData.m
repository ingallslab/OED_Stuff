addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%Section 1: Symbolics for basic functions
x1_sym = SX.sym('x1_sym');
x2_sym = SX.sym('x2_sym');
x_sym=[x1_sym x2_sym];

u1_sym = SX.sym('u1_sym');
u2_sym = SX.sym('u2_sym');
Omega_sym = SX.sym('Omega_sym');

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

g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u2_sym./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u1_sym./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;

g=[g1 g2];
g_x=jacobian(g,x_sym);
g_x_func = Function('g_x_func', {x_sym,u1_sym,u2_sym,theta_sym}, {g_x});

%Section 2: Non-symbolic expressions and parameters
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
Omega=0.01;

theta=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];
par= [theta .5 150];

u0_1=0.25;
u0_2=20;
x1_null=@(x_2,u_2) alpha_1 + beta_1./(1+((x_2./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1);
x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);

g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
g_func2=@(x_2,u_1,u_2) alpha_2 + beta_2*1./(1+(((alpha_1 +beta_1*1./(1+((x_2./K_2)*(1/(1+(u_2./kappa_2).^m_2))).^n_1))./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x_2;

g_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x(1) alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x(2)];

B_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)+x(1) 0;...
    0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)+x(2)];


%Section 3: Data generation
uVals = 0.05:0.1:0.95;
OMEGA=0.01;

phi_1=[0.45 0];
phi_2=[0 45];

Phi=[phi_1;phi_2];


% SSAData=cell(length(uVals),1);
% 
% for i=1:length(uVals)
%     u1 = uVals(i);
%     uV = (phi_2-phi_1)*u1+phi_1;
%     u1 = uV(1);
%     u2 = uV(2);
%     s = dlmread(strcat('2D_Michael/Data/SliceData_Omega=0.01_u=',num2str(u1)),'\t');
%     SSAData{i}=s(1:20,:);
% end
figure;
hold on
NormData=cell(length(uVals),1);
uLow = 18.0019;
uHigh = 20.2516;
for i=1:length(uVals)
    disp(strcat('i = ',num2str(i)));
    
    u1 = uVals(i);
    uV = (phi_2-phi_1)*u1+phi_1;
    u1 = uV(1);
    u2 = uV(2);
    u=norm([u1 u2]);
    g_roots1=@(x_1) g_func1(x_1,u1,u2);
    [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
    
    if u<uLow||u>uHigh
        x2_low=x2_null(x1_low,u1);
        
        A= full(g_x_func([x1_low x2_low],u1,u2,theta));
        B= B_func([x1_low x2_low],u1,u2);
        C = lyap(A,B)/Omega;
        
        NormData{i} = [abs(normrnd(x1_low,sqrt(C(1,1)),[1,20]))' abs(normrnd(x2_low,sqrt(C(2,2)),[1,20]))'];
        scatter3(u1*ones(20,1),u2*ones(20,1),NormData{i}(:,1),'r');
        scatter3(u1*ones(20,1),u2*ones(20,1),NormData{i}(:,2),'b');
    else
        x2_low=x2_null(x1_high,u1);
        x2_high=x2_null(x1_low,u1);
        
        Ah = full(g_x_func([x1_high x2_high],u1,u2,theta));
        Bh = B_func([x1_high x2_low],u1,u2);
        Ch = lyap(Ah,Bh)/Omega;
        
        Al = full(g_x_func([x1_low x2_low],u1,u2,theta));
        Bl = B_func([x1_low x2_high],u1,u2);
        Cl = lyap(Al,Bl)/Omega;
        NormData{i} = [abs(normrnd(x1_low,sqrt(Cl(1,1)),[1,10]))' abs(normrnd(x2_high,sqrt(Cl(2,2)),[1,10]))'];
        NormData{i} = [NormData{i}; [abs(normrnd(x1_high,sqrt(Ch(1,1)),[1,10]))' abs(normrnd(x2_low,sqrt(Ch(2,2)),[1,10]))']];
        scatter3(u1*ones(20,1),u2*ones(20,1),NormData{i}(:,1),'r');
        scatter3(u1*ones(20,1),u2*ones(20,1),NormData{i}(:,2),'b');
    end
    
end
hold off

%Generates symbolics
syms=TwoD_Symbols_MLE(NormData,uVals);

%Generates initial values for constraints
x1_m=[];x2_m=[];
x1_l=[];x2_l=[];
x1_h=[];x2_h=[];

figure;
hold on
for i=1:length(uVals)
    u1 = uVals(i);
    uV = (phi_2-phi_1)*u1+phi_1;
    u1 = uV(1);
    u2 = uV(2);
    u = norm([u1 u2]);
    g_roots1=@(x_1) g_func1(x_1,u1,u2);
    [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000); 
    
    if u<uLow
        x1_m = [x1_m x1_low];
        x2_low=x2_null(x1_high,u1);
        x2_m = [x2_m x2_low];
        scatter3(u1,u2,x1_low,'b');
        scatter3(u1,u2,x2_low,'r');
    elseif u>uHigh
        x1_m = [x1_m x1_high];
        x2_high=x2_null(x1_low,u1);
        x2_m = [x2_m x2_high];
        scatter3(u1,u2, x1_high,'b');
        scatter3(u1,u2, x2_high,'r');
    else
        x1_l=[x1_l x1_low];
        x1_h=[x1_h x1_high];
        x2_low  = x2_null(x1_low,u1);
        x2_high = x2_null(x1_high,u1);
        x2_l=[x2_l x2_low];
        x2_h=[x2_h x2_high];
        scatter3(u1,u2, x1_low,'b');
        scatter3(u1,u2, x1_high,'b');
        scatter3(u1,u2, x2_low,'r');
        scatter3(u1,u2, x2_high,'r');
    end
   
end
hold off
format long g
optimVars = [x1_m'; x1_h'; x1_l'; x2_m'; x2_h'; x2_l'; par']; %These are the variables that are optimized
consts = vertcat(NormData{:}); %These are the fixed parameters in the optimizer
consts = vertcat(consts(:));
consts = [consts;Omega];

llVals=cell(length(par),1);
tarR=cell(length(par),1);
figure;
hold on
obj = @(optV) syms.loglikF(optV,consts);
for j=1:length(par)
    tar=par(j);
    tarR{j} = abs(tar-10:0.05:tar+10);
    for i=1:length(tarR{j})
        Tar=par;
        Tar(j)=tarR{j}(i);
        optV=[x1_m'; x1_h'; x1_l'; x2_m'; x2_h'; x2_l'; Tar'];
        llVals{j}=[llVals{j} obj(optV)];
    end
    subplot(7,2,j);
    plot(tarR{j},llVals{j});
end
hold off
solution = syms.solver('x0',optimVars,'lbg',syms.lbg,'ubg',syms.ubg,'p',consts,'lbx',syms.lbw,'ubx',syms.ubw);
disp(solution);
% disp("hello");
