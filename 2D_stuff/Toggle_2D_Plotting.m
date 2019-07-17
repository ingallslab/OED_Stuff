addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all  

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

u0_1=0.25;
u0_2=20;

% umax_1=0.8;
% umax_2=70;
umax_1=1;
umax_2=100;

x1_null=@(x_2,u_2) alpha_1 + beta_1./(1+((x_2./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1);
x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);

g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
g_func2=@(x_2,u_1,u_2) alpha_2 + beta_2*1./(1+(((alpha_1 +beta_1*1./(1+((x_2./K_2)*(1/(1+(u_2./kappa_2).^m_2))).^n_1))./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x_2;

g_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x(1) alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x(2)];

B_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)+x(1) 0;...
            0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)+x(2)];

%% Casadi Stuff
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

%g_u_func = Function('g_u_func', {x_sym,u_sym,theta_sym}, {g_u});
%g_theta=jacobian(g,theta_sym);
%g1_x=jacobian(g1,x_sym);
%g2_x=jacobian(g2,x_sym);
%g_u=jacobian(g,u_sym);
%full(g_x_func([500 300],u0_1,u0_2,theta))


%% cusp catastrophe plot

[u1_grid,u2_grid] = meshgrid(0:0.01:umax_1,0:1:umax_2);

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


figure
hold on
plot3(xyz_stable_1(:,1),xyz_stable_1(:,2),xyz_stable_1(:,3),'.b')
plot3(xyz_unstable_1(:,1),xyz_unstable_1(:,2),xyz_unstable_1(:,3),'.r')
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


%% 1D input, parameterized line

phi_1=[0.45 0];
phi_2=[0 45];
% phi_1=[0 8];
% phi_2=[0.5 32];
% phi_1=[0.2 0];
% phi_2=[0 20];

x1_vec=[1:4000];
x2_vec=[1:1500];

[u1_grid,u2_grid] = meshgrid(0:0.01:umax_1,0:1:umax_2);

x_tst=[0:1e-2:4000];
Z=zeros(size(u1_grid));
for i=1:size(u1_grid,1)
    i
    for j=1:size(u2_grid,2)
        g_roots1=@(x_1) g_func1(x_1,u1_grid(i,j),u2_grid(i,j));
        gvals=g_roots1(x_tst)+eps;
        sgn_chng=abs(diff(sign(gvals)./2))>0;
        sgn_cnt=sum(sgn_chng);
        if sgn_cnt==1
             Z(i,j)=0;
        else
             Z(i,j)=1;
        end
    end
end
figure
hold on
contourf(u1_grid,u2_grid,Z) 
plot([phi_1(1); phi_2(1)],[phi_1(2); phi_2(2)],'r','LineWidth',3)
hold off
        
%u_vec=[0:0.05:.4 .35:0.01:.55 .6:0.05:1];
%u_vec=[0:0.025:1];
u_vec=[linspace(0,.35,10) linspace(.36,.54,41) linspace(.55,1,10)];
ellps_low={};
ellps_hgh={};
uell_low=[];
uell_hgh=[];

div_surf={};
udiv=[];

STD = 2;                     %# 2 standard deviations
conf = 2*normcdf(STD)-1;     %# covers around 95% of population
scale = chi2inv(conf,2)*5; % ADDED OMEGA FACTOR HERE

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];

x12u=[];
x12u_unstable=[];

u_strt=[.88 59]
%u_strt=[.9 90]
g_roots=@(x_1) g_func1(x_1,u_strt(1),u_strt(2));
x1_strt=fixed_point_v5(g_roots,4000);
x2_strt=x2_null(x1_strt,u_strt(1));
% x2_mid=x2_null(x1_mid,u_vals(1));
% x2_high=x2_null(x1_high,u_vals(1));

A=full(g_x_func([x1_strt x2_strt],u_strt(1),u_strt(2),theta));
B=B_func([x1_strt x2_strt],u_strt(1),u_strt(2));
C = lyap(A,B);
[V D] = eig( scale*C );     %#' cov(X0)
[D order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);
%# unit circle
VV = V*sqrt(D);               %# scale eigenvectors
e_strt = bsxfun(@plus, VV*e, [x1_strt x2_strt]'); %#' project circle back to orig space


figure
N=length(u_vec)-1;
%inc=round(N/16);
cnt=1;
for i=1:N+1
    [i u_vec(i)]
    u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
    g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
    [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
    x2_low=x2_null(x1_low,u_vals(1));
    x2_mid=x2_null(x1_mid,u_vals(1));
    x2_high=x2_null(x1_high,u_vals(1));
    
    x12u=[x12u; [x1_low x2_low u_vec(i)]];

    A_low=full(g_x_func([x1_low x2_low],u_vals(1),u_vals(2),theta));
    B_low=B_func([x1_low x2_low],u_vals(1),u_vals(2));
    C_low = lyap(A_low,B_low);

    [V D] = eig( scale*C_low );     %#' cov(X0)
    [D order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);

    %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e_low = bsxfun(@plus, VV*e, [x1_low x2_low]'); %#' project circle back to orig space

    %# plot cov and major/minor axes
    ellps_low{end+1}=e_low;
    uell_low=[uell_low u_vec(i)];
        
    if ~(x1_low==x1_mid&&x1_mid==x1_high)
        %x12u=[x12u; [x1_low x2_low u_vec(i)]];
        x12u_unstable=[x12u_unstable; [x1_mid x2_mid u_vec(i)]];
        x12u=[x12u; [x1_high x2_high u_vec(i)]];
        
        A_high=full(g_x_func([x1_high x2_high],u_vals(1),u_vals(2),theta));
        B_high=B_func([x1_high x2_high],u_vals(1),u_vals(2));
        C_high = lyap(A_high,B_high);
        
        [V D] = eig( scale*C_high );     %#' cov(X0)
        [D order] = sort(diag(D), 'descend');
        D = diag(D);
        V = V(:, order);

        VV = V*sqrt(D);               %# scale eigenvectors
        e_hgh = bsxfun(@plus, VV*e, [x1_high x2_high]'); %#' project circle back to orig space

        %# plot cov and major/minor axes
        ellps_hgh{end+1}=e_hgh;
        uell_hgh=[uell_hgh u_vec(i)];
        
        J0=full(g_x_func([x1_mid x2_mid],u_vals(1),u_vals(2),theta));
        [eVecs,eVals]=eig(J0);
        eVals=diag(eVals)';
        vectr=eVecs(:,eVals<0);
        vectr=vectr/norm(vectr);
        g_func_fixed = @(t,x) -g_func(x,u_vals(1),u_vals(2))';
        Opt_r = odeset('Events', @myEvent_r);
        Opt_l = odeset('Events', @myEvent_l);
        [~,y_r] = ode45(g_func_fixed,[0 10],[x1_mid; x2_mid]+0.001*vectr,Opt_r);
        [~,y_l] = ode45(g_func_fixed,[0 10],[x1_mid; x2_mid]-0.001*vectr,Opt_l);
        
        div_surf{end+1}=[flipud(y_l); y_r];
        udiv=[udiv u_vec(i)];
        
    end
    if mod(i-1,N/15)==0%(i>=round(mod(N,16)/2)) .* (mod(i-round(mod(N,16)/2),inc)==0) .* (i<(N-(mod(N,16)-round(mod(N,16)/2))))
        subplot(4,4,cnt)
        hold on
        plot(x1_vec,x2_null(x1_vec,u_vals(1)),'-k')
        plot(x1_null(x2_vec,u_vals(2)),x2_vec,'-.k')
        plot(x1_low,x2_low,'b*')
        plot(e_low(1,:), e_low(2,:), 'Color','k');
        if ~(x1_low==x1_mid&&x1_mid==x1_high)
            plot(x1_high,x2_high,'b*')
            plot(x1_mid,x2_mid,'ro')
            quiver(x1_mid,x2_mid,vectr(1),vectr(2),100)
            plot(y_r(:,1),y_r(:,2),'--r')
            plot(y_l(:,1),y_l(:,2),'--r')
            plot([0 4000],(vectr(2)/vectr(1))*([0 4000]-x1_mid)+x2_mid,':b')
            plot(e_hgh(1,:), e_hgh(2,:), 'Color','k');
            plot(e_strt(1,:), e_strt(2,:), 'Color','k');
        end
        xlim([0 4000])
        ylim([0 1500])
        hold off
        cnt=cnt+1;
    end
    
end

figure
hold on
plot3(x12u(:,3),x12u(:,1),x12u(:,2),'ob')
plot3(x12u_unstable(:,3),x12u_unstable(:,1),x12u_unstable(:,2),'.r')
for i=1:length(uell_low)
    e=ellps_low{i};
    plot3(repmat(uell_low(i),1,length(e(1,:))),e(1,:), e(2,:), 'Color','k');
end
for i=1:length(ellps_hgh)
    e=ellps_hgh{i};
    plot3(repmat(uell_hgh(i),1,length(e(1,:))),e(1,:), e(2,:), 'Color','k');
end
for i=1:length(udiv)
    y=div_surf{i};
    plot3(repmat(udiv(i),1,length(y(:,1))),y(:,1)', y(:,2)', 'Color','k');
end
hold off


logLik([500 500],.5,theta,[.5 150],[phi_1;phi_2])


function [value, isterminal, direction] = myEvent_r(T, Y)
    value      = [Y(1)-3500; Y(2)-1250];
    isterminal = [1; 1];   % Stop the integration
    direction  = 0;
end

function [value, isterminal, direction] = myEvent_l(T, Y)
    value      = [Y(1); Y(2)];
    isterminal = [1; 1];   % Stop the integration
    direction  = 0;
end

function ll=logLik(x,u,theta,gamma,Phi)

    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*

    x1_sym = SX.sym('x1_sym');
    x2_sym = SX.sym('x2_sym');
    x_sym=[x1_sym x2_sym];

    u1_sym = SX.sym('u1_sym');
    u2_sym = SX.sym('u2_sym');
    
    u_sym=[u1_sym u2_sym];
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
    g_x_func = Function('g_x_func', {x_sym,u_sym,theta_sym}, {g_x});

    alpha_1 = theta(1);
    alpha_2 =theta(2);
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
    
    phi_1=Phi(1,:);
    phi_2=Phi(2,:);
    u_vals=(phi_2-phi_1)*u+phi_1;

    gamma_0=gamma(1);
    gamma_1=gamma(2);
    
    x1_root=@(x_1) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_vals(1)./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)-x_1;
    x2_eval=@(x_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2);
    B_func=@(x) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1)+x(1) 0;...
            0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2)+x(2)];
    
    [x1_low,~,x1_high]=fixed_point_v5(x1_root,3000);
    x_low=[x1_low x2_eval(x1_low)]
    x_high=[x1_high x2_eval(x1_high)]
    
    A_low=full(g_x_func(x_low,u_vals,theta));
    B_low=B_func(x_low);
    C_low = lyap(A_low,B_low);
    
    A_hgh=full(g_x_func(x_high,u_vals,theta));
    B_hgh=B_func(x_high);
    C_hgh = lyap(A_hgh,B_hgh);
    
    rho=1/(1+exp(-(gamma_1*(u-gamma_0))));
    
    ll=log(rho*mvnpdf(x,x_low,C_low)+(1-rho)*mvnpdf(x,x_high,C_hgh));

end
