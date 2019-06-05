addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%%
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
delta=1.65e-2;

theta=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];


%%

x1_null=@(x_2,u_2) alpha_1 + beta_1./(1+((x_2./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1);
x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);

g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
g_func2=@(x_2,u_1,u_2) alpha_2 + beta_2*1./(1+(((alpha_1 +beta_1*1./(1+((x_2./K_2)*(1/(1+(u_2./kappa_2).^m_2))).^n_1))./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x_2;

g_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x(1) alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x(2)];

B_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)+x(1) 0;...
            0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)+x(2)];

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

%%
OMEGA=.05;
totalTime=20*60;

phi_1=[0.45 0];
phi_2=[0 45];

u=.45;
u_vals=(phi_2-phi_1)*u+phi_1;
%%

% [YOut,tdom]=Toggle_SSA([0 0],u_vals,OMEGA,totalTime);
% 
% figure
% plot(tdom/60,YOut(1,:),tdom/60,YOut(2,:))

x_tst=[0:1e-2:4000];
u_vec=[0:0.01:1];
%Z=zeros(size(u_vec));
u_l=NaN; u_r=NaN;
flag=0;
sgn_cnt_lst=1;
for i=1:length(u_vec)
    u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
    g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
    gvals=g_roots1(x_tst)+eps;
    sgn_chng=abs(diff(sign(gvals)./2))>0;
    sgn_cnt_lst=sgn_cnt;
    sgn_cnt=sum(sgn_chng);
    if (sgn_cnt~=sgn_cnt_lst)&&flag==0
         u_l=u_vec(i);
         flag=1;
    elseif (sgn_cnt~=sgn_cnt_lst)&&flag==1
         u_r=u_vec(i-1);
    end
end

del_u=u_r-u_l;
del_u5p=del_u*.05;

%%

% u_tau=linspace(u_l+del_u5p,u_r-del_u5p,10);
% Omeg_vec=[0.008 0.01 0.02];
% totCnt=100;
% 
% tau_low=zeros(length(u_tau),length(Omeg_vec));
% tau_high=zeros(length(u_tau),length(Omeg_vec));
% 
% for k=1:length(Omeg_vec)
%     Omeg_val=Omeg_vec(k);
%     
%     for i=1:length(u_tau)
%         
%         u_vals=(phi_2-phi_1)*u_tau(i)+phi_1;
% 
%         g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
%         [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
%         x2_low=x2_null(x1_low,u_vals(1));
%         x2_mid=x2_null(x1_mid,u_vals(1));
%         x2_high=x2_null(x1_high,u_vals(1));
%         J0=full(g_x_func([x1_mid x2_mid],u_vals(1),u_vals(2),theta));
%         [eVecs,eVals]=eig(J0);
%         eVals=diag(eVals)';
%         vectr=eVecs(:,eVals<0);
%         vectr=vectr/norm(vectr);
% 
%         vectr=flipud(vectr)';
%         vectr(1)=-1*vectr(1);
% 
%         tau_low_sum=0;
%         tau_high_sum=0;
% 
%         for j=1:totCnt
%             [k i j]
%             tau_low_sum=tau_low_sum+Toggle_SSA_SwitchTime([x1_low x2_low],u_vals,Omeg_val,[x1_mid x2_mid],vectr,1e6);
%             tau_high_sum=tau_high_sum+Toggle_SSA_SwitchTime([x1_high x2_high],u_vals,Omeg_val,[x1_mid x2_mid],vectr,1e6);
%         end
% 
%         tau_low(i,k)=tau_low_sum/totCnt;
%         tau_high(i,k)=tau_high_sum/totCnt;
% 
%     end
% 
% end
% 
% 
% figure
% hold on
% plot(u_tau,tau_low)
% plot(u_tau,tau_high)
% hold off
% 
% 
%   
% %%
%  
% x1_vec=[1:4000];
% x2_vec=[1:1500];
% 
% STD = 2;                     %# 2 standard deviations
% conf = 2*normcdf(STD)-1;     %# covers around 95% of population
% scale = chi2inv(conf,2)*(1/OMEGA); % ADDED OMEGA FACTOR HERE
% 
% t = linspace(0,2*pi,100);
% e = [cos(t) ; sin(t)];
% 
% g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
% [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
% x2_low=x2_null(x1_low,u_vals(1));
% x2_mid=x2_null(x1_mid,u_vals(1));
% x2_high=x2_null(x1_high,u_vals(1));
% 
% A_low=full(g_x_func([x1_low x2_low],u_vals(1),u_vals(2),theta));
% B_low=B_func([x1_low x2_low],u_vals(1),u_vals(2));
% C_low = lyap(A_low,B_low);
% 
% [V D] = eig( scale*C_low );     %#' cov(X0)
% [D order] = sort(diag(D), 'descend');
% D = diag(D);
% V = V(:, order);
% 
% %# unit circle
% VV = V*sqrt(D);               %# scale eigenvectors
% e_low = bsxfun(@plus, VV*e, [x1_low x2_low]'); %#' project circle back to orig space
% 
% if ~(x1_low==x1_mid&&x1_mid==x1_high)
% 
%     A_high=full(g_x_func([x1_high x2_high],u_vals(1),u_vals(2),theta));
%     B_high=B_func([x1_high x2_high],u_vals(1),u_vals(2));
%     C_high = lyap(A_high,B_high);
% 
%     [V D] = eig( scale*C_high );     %#' cov(X0)
%     [D order] = sort(diag(D), 'descend');
%     D = diag(D);
%     V = V(:, order);
% 
%     VV = V*sqrt(D);               %# scale eigenvectors
%     e_hgh = bsxfun(@plus, VV*e, [x1_high x2_high]'); %#' project circle back to orig space
% 
%     J0=full(g_x_func([x1_mid x2_mid],u_vals(1),u_vals(2),theta));
%     [eVecs,eVals]=eig(J0);
%     eVals=diag(eVals)';
%     vectr=eVecs(:,eVals<0);
%     vectr=vectr/norm(vectr);
%     g_func_fixed = @(t,x) -g_func(x,u_vals(1),u_vals(2))';
%     Opt_r = odeset('Events', @myEvent_r);
%     Opt_l = odeset('Events', @myEvent_l);
%     [~,y_r] = ode45(g_func_fixed,[0 10],[x1_mid; x2_mid]+0.001*vectr,Opt_r);
%     [~,y_l] = ode45(g_func_fixed,[0 10],[x1_mid; x2_mid]-0.001*vectr,Opt_l);
% end
% 
% figure
% hold on
% plot(x1_vec,x2_null(x1_vec,u_vals(1)),'-k')
% plot(x1_null(x2_vec,u_vals(2)),x2_vec,'-.k')
% plot(x1_low,x2_low,'b*')
% plot(e_low(1,:), e_low(2,:), 'Color','k');
% if ~(x1_low==x1_mid&&x1_mid==x1_high)
%     plot(x1_high,x2_high,'b*')
%     plot(x1_mid,x2_mid,'ro')
%     quiver(x1_mid,x2_mid,vectr(1),vectr(2),100)
%     plot(y_r(:,1),y_r(:,2),'--r')
%     plot(y_l(:,1),y_l(:,2),'--r')
%     plot([0 4000],(vectr(2)/vectr(1))*([0 4000]-x1_mid)+x2_mid,':b')
%     plot(e_hgh(1,:), e_hgh(2,:), 'Color','k');
% end
% plot(YOut(1,:),YOut(2,:),'g')
% xlim([0 4000])
% ylim([0 1500])
% hold off


%%
num=100;
hrs1=5;
hrs2=20;

u=0;
u_vals=(phi_2-phi_1)*u+phi_1;

g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
[x1_low,~,x1_high]=fixed_point_v5(g_roots1,3000);
x2_low=x2_null(x1_low,u_vals(1));
x2_high=x2_null(x1_high,u_vals(1));

Y_strt=Toggle_Ensemble([x1_low x2_low],u_vals,OMEGA,hrs1*60,num);

u_vec=linspace(u_l+del_u5p,u_r-del_u5p,5);

clrs=lines(num);

figure
for i=1:length(u_vec)
    i
    u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
    
    Yend=Toggle_Ensemble(Y_strt,u_vals,OMEGA,hrs2*60,num);
    
    hold on
    plot(Yend(:,1),Yend(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
    hold off
end
xlim([0 4000])
ylim([0 1500])

u=1;
u_vals=(phi_2-phi_1)*u+phi_1;

g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
[x1_low,~,x1_high]=fixed_point_v5(g_roots1,3000);
x2_low=x2_null(x1_low,u_vals(1));
x2_high=x2_null(x1_high,u_vals(1));

Y_strt=Toggle_Ensemble([x1_low x2_low],u_vals,OMEGA,hrs1*60,num);

figure
for i=1:length(u_vec)
    i
    u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
    
    Yend=Toggle_Ensemble(Y_strt,u_vals,OMEGA,hrs2*60,num);
    
    hold on
    plot(Yend(:,1),Yend(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
    hold off
end
xlim([0 4000])
ylim([0 1500])



%%

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

