addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%FIX SIG, IS WRONG!!!!!

%% Parameter declarations etc.
u_vec=[0:0.005:0.3];
N_u=length(u_vec);
bounds=[0.1 0.2];

%parameters
a0=0.5;
a=3;
K=9;
n=3;
Omega=50; %add this in

theta=[a0 a K n];
N_th=4;

bndryVal=6;
c0=bndryVal*(bounds(1)+bounds(2))/(bounds(1)-bounds(2));
c1=(-bndryVal-c0)/bounds(1);
% c0=-10;
% c1=0;
c=[c0 c1];
N_c=2;

pars=[theta c];
Np=length(pars);

%% Casadi Settup for Liklihood/FIM/Sensitivies

x_sym = SX.sym('x_sym');
u_sym = SX.sym('u_sym');
Omega_sym = SX.sym('Omega_sym');

a0_sym = SX.sym('a0_sym');
a_sym = SX.sym('a_sym');
K_sym = SX.sym('K_sym');
n_sym = SX.sym('n_sym');

theta_sym=[a0_sym a_sym K_sym n_sym];

c0_sym = SX.sym('c0_sym');
c1_sym = SX.sym('c1_sym');

c_sym=[c0_sym c1_sym];

par_sym=[theta_sym c_sym];

g=a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)-x_sym;
%g_func = Function('g_func', {x_sym,u_sym,theta_sym}, {g});

g_theta=jacobian(g,theta_sym);
g_x=jacobian(g,x_sym);
g_u=jacobian(g,u_sym);

g_x_func = Function('g_x_func', {x_sym,u_sym,theta_sym}, {g_x});
g_u_func = Function('g_u_func', {x_sym,u_sym,theta_sym}, {g_u});

%log sensitivities
x_theta=-g_theta./g_x;
%x_theta=(-g_theta./g_x).*theta_sym;

x_theta_func = Function('x_theta_func', {x_sym,u_sym,theta_sym}, {x_theta});

%sigma2=(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)^2/Omega_sym;
sigma2=(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
%log sensitivities
sigma2_theta=jacobian(sigma2,theta_sym);
sigma2_x=jacobian(sigma2,x_sym);
%sigma2_theta=jacobian(sigma2,theta_sym).*theta_sym;

sigma2_func = Function('sigma2_func', {x_sym,u_sym,theta_sym,Omega_sym}, {sigma2});
sigma2_theta_func = Function('sigma2_theta_func', {x_sym,u_sym,theta_sym,Omega_sym}, {sigma2_theta});
sigma2_x_func = Function('sigma2_x_func', {x_sym,u_sym,theta_sym,Omega_sym}, {sigma2_x});

pi0=1/(1+exp(-(c0_sym+c1_sym*u_sym)));
pi0_func = Function('pi0_func', {u_sym,c_sym}, {pi0});

x_obs_sym = SX.sym('x_obs_sym');
phi=(1/(sqrt(2*pi*sigma2)))*exp((-(x_obs_sym-x_sym).^2)./(2*sigma2));

phi_func = Function('phi_func', {x_obs_sym,x_sym,u_sym,theta_sym,Omega_sym}, {phi});

x_low_sym = SX.sym('x_low_sym');
x_high_sym = SX.sym('x_high_sym');

Lik=(1-pi0)*phi_func(x_obs_sym,x_low_sym,u_sym,theta_sym,Omega_sym)+pi0*phi_func(x_obs_sym,x_high_sym,u_sym,theta_sym,Omega_sym);
Lik_func = Function('Lik_func', {x_obs_sym,x_low_sym,x_high_sym,u_sym,par_sym,Omega_sym}, {Lik});
logLik=log(Lik);
logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_sym,par_sym,Omega_sym}, {logLik});

%log sensitivities
logLik_pars=jacobian(logLik,par_sym);
%logLik_pars=jacobian(logLik,par_sym).*par_sym;

logLik_xlow=jacobian(logLik,x_low_sym);
logLik_xhigh=jacobian(logLik,x_high_sym);

logLik_sens=logLik_pars+[logLik_xlow*x_theta_func(x_low_sym,u_sym,theta_sym) zeros(1,N_c)]+[logLik_xhigh.*x_theta_func(x_high_sym,u_sym,theta_sym) zeros(1,N_c)];

logLik_sens_func = Function('logLik_sens_func', {x_obs_sym,x_low_sym,x_high_sym,u_sym,par_sym,Omega_sym}, {logLik_sens});


%% Fitting - Experimental Design Generation

u_bistable=linspace(bounds(1)+.25*(bounds(2)-bounds(1)),bounds(2)-.25*(bounds(2)-bounds(1)),5)';
u_monostable1=linspace(0,bounds(1),3)';
u_monostable2=linspace(bounds(2),.3,3)';

u_opt=[u_monostable1; u_bistable; u_monostable2];
lambda=(1/11)*ones(size(u_opt));


%% Fitting - Plot Bifurcation Diagram

lam_in=lambda;

f1=figure;
f2=figure;

u_mono1=[]; u_mono2=[];u_low=[];u_high=[];
roots_mono1=[]; roots_mono2=[];roots_low=[];roots_high=[];
flg=0;
for i=1:length(u_vec)
    [lPt,mPt,hPt]=fixed_point_v3(u_vec(i),theta);

    if (lPt==hPt&&flg==0)
        u_mono1=[u_mono1 u_vec(i)];
        roots_mono1=[roots_mono1 lPt];

    elseif(lPt==hPt&&flg==1)
        u_mono2=[u_mono2 u_vec(i)];
        roots_mono2=[roots_mono2 lPt];
    else
        u_low=[u_low u_vec(i)];
        roots_low=[roots_low lPt];
        u_high=[u_high u_vec(i)];
        roots_high=[roots_high hPt];
        flg=1;
    end
end

figure(f1);
hold on
plot(u_mono1,roots_mono1,'Color',[0 0 0]);
plot(u_mono2,roots_mono2,'Color',[0 0 0]);
plot(u_low,roots_low,'Color',[0 0 0]);
plot(u_high,roots_high,'Color',[0 0 0]);
hold off

ys=1./(1+exp(-(c0+c1*u_vec)));
figure(f2);
hold on
plot(u_vec,ys,'Color',[0 0 0]);
hold off

%% Fitting - Data Generation

lam_round=lam_in(lam_in>1e-4);
%u_opt=u_vec(lam_in>1e-4);
num_samps=round(lam_round*1000);

Data=[];
Input=[];

DataCell={};

for i=1:length(u_opt)
    [lPt,mPt,hPt]=fixed_point_v3(u_opt(i),theta);
    pi_0=full(pi0_func(u_opt(i),c));
    up_N=round(num_samps(i)*pi_0);
    dwn_N=num_samps(i)-up_N;
    
    %dwn_samp=normrnd(lPt,sqrt(full(sigma2_func(lPt,u_opt(i),theta,Omega))),dwn_N,1);
    %up_samp=normrnd(hPt,sqrt(full(sigma2_func(hPt,u_opt(i),theta,Omega))),up_N,1);
    %dwn_samp=normrnd(lPt,.7*sqrt(full(sigma2_func(lPt,u_opt(i),theta,Omega))),dwn_N,1);
    %up_samp=normrnd(hPt,.7*sqrt(full(sigma2_func(hPt,u_opt(i),theta,Omega))),up_N,1);
    dwn_samp=repmat(lPt,dwn_N,1);%normrnd(lPt,sqrt(full(sigma2_func(lPt,u_opt(i),theta,Omega))),dwn_N,1);
    up_samp=repmat(hPt,up_N,1);%normrnd(hPt,sqrt(full(sigma2_func(hPt,u_opt(i),theta,Omega))),up_N,1);

    samps=[dwn_samp;up_samp];
    samps(samps<0)=0;
    
    Data=[Data; samps];
    Input=[Input; repmat(u_opt(i),num_samps(i),1)];
    
    DataCell{i}=samps;
end

N_obs=sum(num_samps);
Tvec=zeros(N_obs,1);
Tcell={};

theta_0=normrnd(theta,abs(.1*theta));
%theta_0=normrnd(theta,abs(.25*theta));
%theta_0=[.4 2.8 8 3.2];
%theta_0=[.4 3.1 9 3.2];
%theta_0=[.3 2.2 4.4 3.7];

par_est=[theta_0 c];
par_est_last=zeros(size(pars));

clr=spring(3);
cntr=1;

figure(f1);
hold on
plot(Input,Data,'*r');
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the EM Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% while norm(par_est-par_est_last)>.01
%     
%     theta_est=par_est(1:N_th);
%     c_est=par_est(N_th+1:end);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Plotting, turn off later
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     u_mono1=[]; u_mono2=[];u_low=[];u_high=[];
%     roots_mono1=[]; roots_mono2=[];roots_low=[];roots_high=[];
%     flg=0;
%     for i=1:length(u_vec)
%         [lPt,mPt,hPt]=fixed_point_v3(u_vec(i),par_est(1:N_th));
% 
%         if (lPt==hPt&&flg==0)
%             u_mono1=[u_mono1 u_vec(i)];
%             roots_mono1=[roots_mono1 lPt];
%             
%         elseif(lPt==hPt&&flg==1)
%             u_mono2=[u_mono2 u_vec(i)];
%             roots_mono2=[roots_mono2 lPt];
%         else
%             u_low=[u_low u_vec(i)];
%             roots_low=[roots_low lPt];
%             u_high=[u_high u_vec(i)];
%             roots_high=[roots_high hPt];
%             flg=1;
%         end
%     end
%     figure(f1)
%     hold on
%     plot(u_mono1,roots_mono1,u_mono2,roots_mono2,...
%         u_low,roots_low,u_high,roots_high,'Color',clr(cntr,:));
%     hold off
% 
%     ys=1./(1+exp(-(c_est(1)+c_est(2)*u_vec)));
%     figure(f2);
%     hold on
%     plot(u_vec,ys,'Color',clr(cntr,:));
%     hold off
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expectation Step
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cnt1=1;
%     for i=1:length(u_opt)
%         N_uopt=length(DataCell{i});
%         
%         [lPt,mPt,hPt]=fixed_point_v3(u_opt(i),theta_est);
%         std_dwn=sqrt(full(sigma2_func(lPt,u_opt(i),theta_est,Omega)));
%         std_up=sqrt(full(sigma2_func(hPt,u_opt(i),theta_est,Omega)));
%         
%         if(u_opt(i)<=bounds(1))
%             pi_0=0;
%         elseif(u_opt(i)>=bounds(2))
%             pi_0=1;
%         else
%             pi_0=1./(1+exp(-(c_est(1)+c_est(2).*u_opt(i))));%full(pi0_func(u_opt(i),c_est));
%         end
%         
%         Tvec(cnt1:cnt1+(N_uopt-1))=(pi_0.*normpdf(DataCell{i},hPt,std_up))...
%                                     ./((1-pi_0).*normpdf(DataCell{i},lPt,std_dwn)...
%                                         +pi_0.*normpdf(DataCell{i},hPt,std_up));
%         Tcell{i}=Tvec(cnt1:cnt1+(N_uopt-1));
%         cnt1=cnt1+N_uopt;
%     end
%     
%     par_est_last=par_est;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Maximization Step
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Fit logistic, c pars
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %this is a hack for now, just to use matlab's logistic fitting alg.
%     %which requires counts (hence 100*Tvec) rather than proportions (just
%     %Tvec)
%     c_est=glmfit(Input,[round(Tvec*100) 100*ones(size(Tvec))],'binomial','logit');
%     
%     % Fit steady state model, theta pars
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fit_func=@(th) fit_obj(th,u_opt,DataCell,Tcell,bounds,Omega);
%     options = optimoptions('particleswarm','Display','iter','HybridFcn',@fmincon,'FunctionTolerance',1e-3);
%     fit_func(theta)
%     theta_est = particleswarm(fit_func,N_th,[0.001 0.001 1 1],[1 5 10 4],options);
%     
%     par_est=[theta_est c_est'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cntr=cntr+1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try One-step Max Liklihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit_func=@(prs) logLikObj(prs,u_opt,DataCell,bounds,Omega);
options = optimoptions('particleswarm','Display','iter','SwarmSize',200,'FunctionTolerance',1e-4);
fit_func(pars)
pars_est = particleswarm(fit_func,Np,[0.001 0.001 1 1 -Inf -Inf],[1 5 10 4 Inf Inf],options);


%% Functions

function err=fit_obj(theta,u_opt,DataCell,Tcell,Bounds,Omega)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    tol=1e-4;
    %x_tst=[0:tol:maxU];

    g_bnd1_cnt=0;
    maxU=5;
    while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
        x_tst=[0:tol:maxU];
        g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    
    g_bnd2_cnt=0;
    maxU=5;
    while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
        x_tst=[0:tol:maxU];
        g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end

%     g_bnd1_cnt=sum(abs(diff(sign(a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)./2)));
%     g_bnd2_cnt=sum(abs(diff(sign(a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)./2)));
      
    if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
        err=1e8;
        return
    end
    
    err=0;
    for i=1:length(u_opt)
        [lPt,~,hPt]=fixed_point_v3(u_opt(i),theta);

        if(u_opt(i)<=Bounds(1))
            sig_low2=((a0+a*((u_opt(i)+lPt)^n)/(K+(u_opt(i)+lPt)^n)+lPt).^2/Omega);
            err=err+sum(((1-Tcell{i}).*((DataCell{i}-lPt)).^2)./sig_low2);
        elseif(u_opt(i)>=Bounds(2))
            sig_hgh2=((a0+a*((u_opt(i)+hPt)^n)/(K+(u_opt(i)+hPt)^n)+hPt).^2/Omega);
            err=err+sum((Tcell{i}).*((DataCell{i}-hPt).^2)./sig_hgh2);
        else
            sig_low2=((a0+a*((u_opt(i)+lPt)^n)/(K+(u_opt(i)+lPt)^n)+lPt).^2/Omega);
            sig_hgh2=((a0+a*((u_opt(i)+hPt)^n)/(K+(u_opt(i)+hPt)^n)+hPt).^2/Omega);
            err=err+sum((1-Tcell{i}).*((DataCell{i}-lPt).^2)./sig_low2+(Tcell{i}).*((DataCell{i}-hPt).^2)./sig_hgh2);
        end
    end
    
end

function [low,mid,high]=fixed_point_v3(u,theta)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    tol=1e-4;
    maxU=5;
    %x_tst=[0:tol:maxU];
    
    sgn_cnt=0;
    while sgn_cnt~=1&&sgn_cnt~=3
        x_tst=[0:tol:maxU];
        gvals=(a0+a.*((u+x_tst).^n)./(K+(u+x_tst).^n)-x_tst)+eps;
        sgn_chng=abs(diff(sign(gvals)./2))>0;
        sgn_cnt=sum(sgn_chng);
        maxU=maxU+1;
    end
    
    x_strt_pts=(x_tst([sgn_chng false]).*gvals([false sgn_chng])...
                -x_tst([false sgn_chng]).*gvals([sgn_chng false]))...
                    ./(gvals([false sgn_chng])-gvals([sgn_chng false]));
    %x_strt_pts=unique(x_strt_pts);
    
    if ~(sgn_cnt==1||sgn_cnt==3)
        low=NaN;
        mid=NaN;
        high=NaN;
        return 
    end

    func=@(x) g_func(x,u,theta);
    low=fzero(func,x_strt_pts(1));
    if(sgn_cnt==1)
        mid=low;
        high=low;
    else
        mid=fzero(func,x_strt_pts(2));
        high=fzero(func,x_strt_pts(3));
    end
end


function ll=logLikObj(pars,u_opt,DataCell,Bounds,Omega)

    theta=pars(1:4);

    a0=pars(1);
    a=pars(2);
    K=pars(3);
    n=pars(4);
    
    c_0=pars(5);
    c_1=pars(6);

    tol=1e-4;
    %x_tst=[0:tol:maxU];

    g_bnd1_cnt=0;
    maxU=5;
    while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
        x_tst=[0:tol:maxU];
        g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    
    g_bnd2_cnt=0;
    maxU=5;
    while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
        x_tst=[0:tol:maxU];
        g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end

%     g_bnd1_cnt=sum(abs(diff(sign(a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)./2)));
%     g_bnd2_cnt=sum(abs(diff(sign(a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)./2)));
      
    if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
        ll=1e8;
        return
    end
    
    ll=0;
    for i=1:length(u_opt)
        [lPt,~,hPt]=fixed_point_v3(u_opt(i),theta);

        if(u_opt(i)<=Bounds(1))
            sig_low=(a0+a*((u_opt(i)+lPt)^n)/(K+(u_opt(i)+lPt)^n)+lPt)/sqrt(Omega);
            %ll=ll+sum(-((DataCell{i}-lPt).^2)./(2*sig_low))-(0.5)*log(2*pi*sig_low);
            ll=ll+sum(log(normpdf(DataCell{i},lPt,sig_low)));
        elseif(u_opt(i)>=Bounds(2))
            sig_hgh=(a0+a*((u_opt(i)+hPt)^n)/(K+(u_opt(i)+hPt)^n)+hPt)/sqrt(Omega);
            %ll=ll+sum(-((DataCell{i}-hPt).^2)./(2*sig_hgh))-(0.5)*log(2*pi*sig_low);
            ll=ll+sum(log(normpdf(DataCell{i},hPt,sig_hgh)));
        else
            sig_low=(a0+a*((u_opt(i)+lPt)^n)/(K+(u_opt(i)+lPt)^n)+lPt)/sqrt(Omega);
            sig_hgh=(a0+a*((u_opt(i)+hPt)^n)/(K+(u_opt(i)+hPt)^n)+hPt)/sqrt(Omega);
            pi_0=1./(1+exp(-(c_0+c_1.*u_opt(i))));
            %ll=ll+sum((1-pi_0).*exp(-((DataCell{i}-lPt).^2)./(2*sig_low))+(pi_0).*exp(-((DataCell{i}-hPt).^2)./sig_hgh));
            ll=ll+sum(log((1-pi_0).*normpdf(DataCell{i},lPt,sig_low)+pi_0.*normpdf(DataCell{i},hPt,sig_hgh)));
        end
    end
    
    ll=-1*ll;
end


function val=g_func(x0,u,theta)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    val=a0+a*((u+x0).^n)./(K+(u+x0).^n)-x0;

end

% function sig=sig_eval(x,u,theta,Omega)
% 
%     a0=theta(1);
%     a=theta(2);
%     K=theta(3);
%     n=theta(4);
% 
%     sig=(a0+a*((u+x)^n)/(K+(u+x)^n)+x)/sqrt(Omega);
% 
% end

