umax_1=1;
umax_2=100;
uVals = 0.05:0.05:0.95;
OMEGA=0.1;

phi_1=[0.45 0];
phi_2=[0 45];

Phi=[phi_1;phi_2];

finTime = 5000;
num = 500;


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
[u1_grid,u2_grid] = meshgrid(0:0.05:umax_1,0:2:umax_2);

xyz_stable_1=[];
xyz_stable_2=[];
xyz_unstable_1=[];
xyz_unstable_2=[];

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
        else
            xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_low]];
            xyz_unstable_1=[xyz_unstable_1; [u1_grid(i,j) u2_grid(i,j) x1_mid]];
            xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_high]];
            
            x2_low=x2_null(x1_low,u1_grid(i,j));
            x2_mid=x2_null(x1_mid,u1_grid(i,j));
            x2_high=x2_null(x1_high,u1_grid(i,j));
            
            xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_low]];
            xyz_unstable_2=[xyz_unstable_2; [u1_grid(i,j) u2_grid(i,j) x2_mid]];
            xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_high]];
            
        end
        
    end
end

xyz = [];
for i=1:length(uVals)
    u1 = uVals(i);
    uV = (phi_2-phi_1)*u1+phi_1;
    u1 = uV(1);
    u2 = uV(2);
    s = dlmread(strcat('2D_Michael/Data/SliceData_Omega=0.01_u=',num2str(u1)),'\t');
    xyz = [xyz; [u1*ones(1,5000)' u2*ones(1,5000)' s(:,1) s(:,2)]];
end

x1_m=[];x2_m=[];
x1_l=[];x2_l=[];
x1_h=[];x2_h=[];
uLow = 18.0019;
uHigh = 20.2516;
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
        x2_low=x2_null(x1_low,u1);
        x2_m = [x2_m x2_low];
        scatter3(u1,u2,x1_low);
    elseif u>uHigh
        x1_m = [x1_m x1_high];
        x2_high=x2_null(x1_high,u1);
        x2_m = [x2_m x2_high];
        scatter3(u1,u2, x1_high);
    else
        x1_l=[x1_l x1_low];
        x1_h=[x1_h x1_high];
        x2_high=x2_null(x1_low,u1);
        x2_low=x2_null(x1_high,u1);
        x2_l=[x2_l x2_low];
        x2_h=[x2_h x2_high];
        scatter3(u1,u2, x1_low);
        scatter3(u1,u2, x1_high);
    end
    
end
plot3(xyz_stable_1(:,1),xyz_stable_1(:,2),xyz_stable_1(:,3),'.b')
plot3(xyz_unstable_1(:,1),xyz_unstable_1(:,2),xyz_unstable_1(:,3),'.b')
hold off

figure
hold on
plot3(xyz_stable_1(:,1),xyz_stable_1(:,2),xyz_stable_1(:,3),'.b')
plot3(xyz_unstable_1(:,1),xyz_unstable_1(:,2),xyz_unstable_1(:,3),'.b')
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'.b')
xlim([0 umax_1])
ylim([0 umax_2])
hold off
figure
hold on
plot3(xyz_stable_2(:,1),xyz_stable_2(:,2),xyz_stable_2(:,3),'.r')
plot3(xyz_unstable_2(:,1),xyz_unstable_2(:,2),xyz_unstable_2(:,3),'.r')
scatter3(xyz(:,1),xyz(:,2),xyz(:,4),'.r')
xlim([0 umax_1])
ylim([0 umax_2])
hold off
