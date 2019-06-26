umax_1=1;
umax_2=100;
[u1_grid,u2_grid] = meshgrid(0:0.05:umax_1,0:5:umax_2);
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
numVals=100;

addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all
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



xyz_stable_1=[];
xyz_stable_2=[];
xyz_unstable_1=[];
xyz_unstable_2=[];
XYZ_stable = [];
XYZ_unstable_high = [];
XYZ_unstable_low = [];
for i=1:size(u1_grid,1)
    disp(strcat('i = ',num2str(i)));
    for j=1:size(u2_grid,2)
        %[i j]
        g_roots1=@(x_1) g_func1(x_1,u1_grid(i,j),u2_grid(i,j));
        [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
        
        if (x1_low==x1_mid&&x1_mid==x1_high)
            xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_low]];
            
            x2_low=x2_null(x1_low,u1_grid(i,j));
            xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_low]];
            
            A=full(g_x_func([x1_low x2_low],u1_grid(i,j),u2_grid(i,j),theta));
            B=B_func([x1_low x2_low],u1_grid(i,j),u2_grid(i,j));
            C = lyap(A,B);
            
            XYZ_stable = [XYZ_stable; [u1_grid(i,j)*ones(1,numVals)' u2_grid(i,j)*ones(1,numVals)' abs(normrnd(x1_low,sqrt(C(1,1)),[1,numVals]))']];
            XYZ_stable = [XYZ_stable; [u1_grid(i,j)*ones(1,numVals)' u2_grid(i,j)*ones(1,numVals)' abs(normrnd(x1_high,sqrt(C(1,1)),[1,numVals]))']];
        else
            xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_low]];
            xyz_unstable_1=[xyz_unstable_1; [u1_grid(i,j) u2_grid(i,j) x1_mid]];
            xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_high]];
            
            x2_low=x2_null(x1_low,u1_grid(i,j));
            x2_mid=x2_null(x1_mid,u1_grid(i,j));
            x2_high=x2_null(x1_high,u1_grid(i,j));
            
            A=full(g_x_func([x1_low x2_low],u1_grid(i,j),u2_grid(i,j),theta));
            B=B_func([x1_low x2_low],u1_grid(i,j),u2_grid(i,j));
            C = lyap(A,B);
            
            xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_low]];
            xyz_unstable_2=[xyz_unstable_2; [u1_grid(i,j) u2_grid(i,j) x2_mid]];
            xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_high]];
            XYZ_stable = [XYZ_stable; [u1_grid(i,j)*ones(1,numVals)' u2_grid(i,j)*ones(1,numVals)' abs(normrnd(x1_low,sqrt(C(1,1)),[1,numVals]))']];
            XYZ_stable = [XYZ_stable; [u1_grid(i,j)*ones(1,numVals)' u2_grid(i,j)*ones(1,numVals)' abs(normrnd(x1_high,sqrt(C(1,1)),[1,numVals]))']];
        end
        
    end
end

figure
hold on
scatter3(XYZ_stable(:,1),XYZ_stable(:,2),XYZ_stable(:,3),'.b')

hold off
xlim([0 umax_1])
ylim([0 umax_2])
figure
hold on
plot3(xyz_stable_2(:,1),xyz_stable_2(:,2),xyz_stable_2(:,3),'.b')
plot3(xyz_unstable_2(:,1),xyz_unstable_2(:,2),xyz_unstable_2(:,3),'.r')
hold off
xlim([0 umax_1])
ylim([0 umax_2])
