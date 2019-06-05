addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all
%figure(1); clf; hold on; box on;

Nh = 2000;
Nensmb = 100000;
%Nensmb = 600000;

%parameters
a0=0.5;
a=3;
K=9;
n=3;
Omega=90; %add this in

theta=[a0 a K n];

u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)] ;%linspace(0.2,0.3,32)
totU=length(u);

xmin = 0;
xmax = 4.0;
histgrid = linspace(xmin,xmax,Nh);

xhist=[];
uTot_new=0;
u_new=[];
for i=1:totU
    [i totU]
    %strcat('drive_2/hist_W=90.000000_u=',num2str(u(i),'%1.6f'),'.txt')
    %tmp=load(strcat('data/hist_W=90.000000_u=',num2str(u(i),'%1.6f'),'_Nr=100000.txt'));
    tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u(i),'%1.6f'),'.txt'));
    sz=size(tmp);
    if i==1||sz(1)==sz_hst(1)
        xhist = [xhist tmp];
        sz_hst = size(xhist);
        uTot_new = uTot_new+1;
        u_new=[u_new u(i)];
    else
        testy=0;
    end
    
end

totU=uTot_new;
u=u_new;


u_W120 = [ linspace(0.12,0.19,12) ] ;
u_W120=[u_W120(1:6) u_W120(8:12) ];
totU_W120=length(u_W120);

xhist_W120={};
for i=1:totU_W120
    [i totU_W120]
    xhist_W120{i} = load(strcat('drive_W120/hist_W=120.000000_u=',num2str(u_W120(i),'%1.6f'),'_Nr=100000.txt'));
end

u_W60 = [ linspace(0.1,0.2,36) ] ;%linspace(0.2,0.3,32)
totU_W60=length(u_W60);

xhist_W60={};
for i=1:totU_W60
    [i totU_W60]
    xhist_W60{i} = load(strcat('drive_W60/hist_W=60.000000_u=',num2str(u_W60(i),'%1.6f'),'.txt'));
end

uL=continAlg(theta,0.1,-1,2)
uR=continAlg(theta,0.2,1,1)

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

sigma2=-(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
%log sensitivities
sigma2_theta=jacobian(sigma2,theta_sym);
sigma2_x=jacobian(sigma2,x_sym);
%sigma2_theta=jacobian(sigma2,theta_sym).*theta_sym;

sigma2_func = Function('sigma2_func', {x_sym,u_sym,theta_sym,Omega_sym}, {sigma2});


sig_EK=(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(Omega_sym);
sigma_EK_func = Function('sigma_EK_func', {x_sym,u_sym,theta_sym,Omega_sym}, {sig_EK});

V_d=a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)-x_sym;
V_dd_sym=-jacobian(V_d,x_sym);

V_dd_func= Function('V_dd', {x_sym,u_sym,theta_sym}, {V_dd_sym});

%addisons code block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kern_pdf={};
kern_inpt={};
for i=1:totU
    %xhistc(:,i) = histc(xhist(:,i),histgrid)/sum(xhist(:,i));ksdensity(x)
    [df,xi]=ksdensity(xhist(:,i));
    kern_pdf{i}=df;
    kern_inpt{i}=xi;
end

u_vec=[0:0.005:0.325];
u_mono1=[]; u_mono2=[];u_low=[];u_high=[];
roots_mono1=[]; roots_mono2=[];roots_low=[];roots_high=[];
sig_mono1=[]; sig_mono2=[]; sig_low=[]; sig_high=[];
u_mid=[];roots_mid=[];
flg=0;
for i=1:length(u_vec)
    [lPt,mPt,hPt]=fixed_point_v3(u_vec(i),theta);

    if (lPt==hPt&&flg==0)
        u_mono1=[u_mono1 u_vec(i)];
        roots_mono1=[roots_mono1 lPt];
        %sig_mono1=[sig_mono1 ((a0+a*((u_vec(i)+lPt)^n)/(K+(u_vec(i)+lPt)^n)+lPt)/sqrt(Omega))];
        sig_mono1=[sig_mono1 sqrt(full(sigma2_func(lPt,u_vec(i),theta,Omega)))];

    elseif(lPt==hPt&&flg==1)
        u_mono2=[u_mono2 u_vec(i)];
        roots_mono2=[roots_mono2 lPt];
        %sig_mono2=[sig_mono2 ((a0+a*((u_vec(i)+lPt)^n)/(K+(u_vec(i)+lPt)^n)+lPt)/sqrt(Omega))];
        sig_mono2=[sig_mono2 sqrt(full(sigma2_func(lPt,u_vec(i),theta,Omega)))];
    else
        u_mid=[u_mid u_vec(i)];
        roots_mid=[roots_mid mPt];
        
        u_low=[u_low u_vec(i)];
        roots_low=[roots_low lPt];
        %sig_low=[sig_low ((a0+a*((u_vec(i)+lPt)^n)/(K+(u_vec(i)+lPt)^n)+lPt)/sqrt(Omega))];
        sig_low=[sig_low sqrt(full(sigma2_func(lPt,u_vec(i),theta,Omega)))];
        u_high=[u_high u_vec(i)];
        roots_high=[roots_high hPt];
        %sig_high=[sig_high ((a0+a*((u_vec(i)+hPt)^n)/(K+(u_vec(i)+hPt)^n)+hPt)/sqrt(Omega))];
        sig_high=[sig_high sqrt(full(sigma2_func(hPt,u_vec(i),theta,Omega)))];
        flg=1;
    end
end

U_Low=[u_mono1 u_low];
X_Low=[roots_mono1 roots_low];
Sig_Low=[sig_mono1 sig_low];
U_Hgh=[u_high u_mono2];
X_Hgh=[roots_high roots_mono2];
Sig_Hgh=[sig_high sig_mono2];

f1=figure; %clf; 
%set(gcf, 'InvertHardcopy', 'off')
subplot(7,1,4:7)
text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
hold on
plot(U_Low,X_Low,'b');
plot(u_mid,roots_mid,'k:');
plot(U_Hgh,X_Hgh,'r');
h1=fill([U_Low';flipud(U_Low')],[X_Low'-Sig_Low';flipud(X_Low'+Sig_Low')],[.4 .4 .7],'linestyle','none');%,'-','LineWidth',2);
set(h1,'facealpha',.3)
h2=fill([U_Hgh';flipud(U_Hgh')],[X_Hgh'-Sig_Hgh';flipud(X_Hgh'+Sig_Hgh')],[.7 .4 .4],'linestyle','none');%,'-','LineWidth',2);
set(h2,'facealpha',.3)
hold off
%xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Concentration','FontSize',14,'FontWeight','bold','Color','k')
%title('Bifurcation Plots for Sample of \theta','FontSize',20,'FontWeight','bold','Color','k')
xlim([0 0.325]);

gap_1=10;

u_1=u(u<0.07);
xhist_1=xhist(:,u<0.07);
uNum_1=length(u_1);

MeanVec=[];
StdVec=[];
for i=5:gap_1:uNum_1 
    [~,mPt,~]=fixed_point_v3(u_1(i),theta);
    MeanVec=[MeanVec mean(xhist_1(:,i))];
    StdVec=[StdVec std(xhist_1(:,i))];
end
hold on
errorbar(u_1(5:gap_1:uNum_1),MeanVec,StdVec,'o','Color',[.2 .2 .3 ]);
%errorbar(u(1:gap:totU),hghMeanVec,hghStdVec,'--','Color',[.3 .2 .2 ]);
hold off

gap_2=7;

u_2=u(u>0.22);
xhist_2=xhist(:,u>0.22);
uNum_2=length(u_2);

MeanVec=[];
StdVec=[];
for i=3:gap_2:uNum_2 
    [~,mPt,~]=fixed_point_v3(u_2(i),theta);
    MeanVec=[MeanVec mean(xhist_2(:,i))];
    StdVec=[StdVec std(xhist_2(:,i))];
end
hold on
errorbar(u_2(3:gap_2:uNum_2),MeanVec,StdVec,'o','Color',[.2 .2 .3 ]);
%errorbar(u(1:gap:totU),hghMeanVec,hghStdVec,'--','Color',[.3 .2 .2 ]);
hold off

gap_3=20;

u_3=u(and(0.1<u,u<0.2));
xhist_3=xhist(:,and(0.1<u,u<0.2));
uNum_3=length(u_3);

lowMeanVec=[];
hghMeanVec=[];
lowStdVec=[];
hghStdVec=[];
for i=1:gap_3:uNum_3 
    [~,mPt,~]=fixed_point_v3(u_3(i),theta);
    lowMeanVec=[lowMeanVec mean(xhist_3(xhist_3(:,i)<mPt,i))];
    hghMeanVec=[hghMeanVec mean(xhist_3(xhist_3(:,i)>mPt,i))];
    lowStdVec=[lowStdVec std(xhist_3(xhist_3(:,i)<mPt,i))];
    hghStdVec=[hghStdVec std(xhist_3(xhist_3(:,i)>mPt,i))];
end
hold on
errorbar(u_3(1:gap_3:uNum_3),lowMeanVec,lowStdVec,'o','Color',[.2 .2 .3 ]);
errorbar(u_3(1:gap_3:uNum_3),hghMeanVec,hghStdVec,'o','Color',[.3 .2 .2 ]);
hold off

legend('Lower Branch','Unstable Points','Upper Branch','LNA Lower Branch','LNA Upper Branch','SSA Mean & StD','Location', 'southeast')

subplot(7,1,1:2)
text(-.16,1.4,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
hold on; %box on;
gap=20;
phnls=[];
for i=1:gap:totU 
    phnls=[phnls plot(kern_inpt{i},kern_pdf{i}, '-','Color', [i/(totU), 0, (totU-i)/(totU)], 'LineWidth', 2.0)];
end
hold off

%legend([phnls(1) phnls(end)],'u=0.0','u=0.3','Location', 'northeast');
colorMap=[(1:totU)'/(totU), zeros(totU,1), (totU-(1:totU)')/(totU)];
colormap(colorMap);
h=colorbar('east');
caxis([0 0.3]);
ylabel(h, 'Input, u')
xlim([0 4])
ylim([0 4.5])
xlabel('Simulated Concentration, X/\Omega','FontSize',14,'FontWeight','bold','Color','k');
ylabel('Estimated PDF','FontSize',14,'FontWeight','bold','Color','k');
%title('PDF in Bistable Region');

saveas(f1,'Scurve_Hist.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some of my stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ntot=length(xhist(:,1));
%Nsamp=10000;

SampVec=[100 1000 10000];
NumSampVals=length(SampVec);

f2=figure
%set(gcf, 'InvertHardcopy', 'off')
subplot(2,1,2)
text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
for j=1:NumSampVals

    HS_Rslt=zeros(NumSampVals,totU);
    %BC_Rslt=zeros(1,12);
    for i=1:totU
       subset=randi(SampVec(j),Ntot,1);
       %BC_Rslt(i)=bmtest(xhist(subset,i));
       [HS_Rslt(j,i)]=HartigansDipTest(sort(xhist(subset,i)));%HartigansDipSignifTest(sort(xhist(subset,i)),500);%;
    end

    hold on
    plot(u,HS_Rslt(j,:))
    hold off
end
xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
ylabel('HDS','FontSize',14,'FontWeight','bold','Color','k')
labs=cellstr(num2str(SampVec', 'N=%-d'));
legend(labs,'Location','northwest');
xlim([0.1 .22])



u_3=u(and(0.075<u,u<0.21));
xhist_3=xhist(:,and(0.075<u,u<0.21));
%uNum_3=length(u_3);

u_3=u_3(1:5:end);
xhist_3=xhist_3(:,1:5:end);
uNum_3=length(u_3);

pi0=zeros(size(uNum_3));
piVal=zeros(size(uNum_3 ));
for i=1:uNum_3 
    [lPt,mPt,hPt]=fixed_point_v3(u_3(i),theta);
    pi0(i) = sum(xhist_3(:,i)>mPt)/sum(xhist_3(:,i)>=0);
    
    
    func=@(x) -g_func(x,u_3(i),theta);
    delt_V_up=integral(func,hPt,mPt);
    delt_V_dwn=integral(func,lPt,mPt);

%     sig2_up=full(sigma2_func(hPt,u_3(i),theta,Omega));
%     sig2_dwn=full(sigma2_func(lPt,u_3(i),theta,Omega));
    
    sig2_up=full(sigma_EK_func(hPt,u_3(i),theta,Omega));
    sig2_dwn=full(sigma_EK_func(lPt,u_3(i),theta,Omega));

    V_dd_up=full(V_dd_func(hPt,u_3(i),theta));
    V_dd_dwn=full(V_dd_func(lPt,u_3(i),theta));

    expVal=2*(delt_V_dwn/sig2_dwn- delt_V_up/sig2_up)+(1/2)*(log(V_dd_up)-log(V_dd_dwn));
    %expVal=2*(delt_V_dwn/sig2_dwn - delt_V_up/sig2_up)+(1/2)*(log(V_dd_up_tst)-log(V_dd_dwn_tst));
    piVal(i)=1/(1+exp(expVal));
end


pi0_W60=zeros(size(totU_W60));
piVal_W60=zeros(size(totU_W60 ));
for i=1:totU_W60 
    [lPt,mPt,hPt]=fixed_point_v3(u_W60(i),theta);
    pi0_W60(i) = sum(xhist_W60{i}>mPt)/sum(xhist_W60{i}>=0);
    
    
    func=@(x) -g_func(x,u_W60(i),theta);
    delt_V_up=integral(func,hPt,mPt);
    delt_V_dwn=integral(func,lPt,mPt);

%     sig2_up=full(sigma2_func(hPt,u_W60(i),theta,60));
%     sig2_dwn=full(sigma2_func(lPt,u_W60(i),theta,60));
    
    sig2_up=full(sigma_EK_func(hPt,u_W60(i),theta,60));
    sig2_dwn=full(sigma_EK_func(lPt,u_W60(i),theta,60));

    V_dd_up=full(V_dd_func(hPt,u_W60(i),theta));
    V_dd_dwn=full(V_dd_func(lPt,u_W60(i),theta));

    expVal=2*(delt_V_dwn/sig2_dwn- delt_V_up/sig2_up)+(1/2)*(log(V_dd_up)-log(V_dd_dwn));
    %expVal=2*(delt_V_dwn/sig2_dwn - delt_V_up/sig2_up)+(1/2)*(log(V_dd_up_tst)-log(V_dd_dwn_tst));
    piVal_W60(i)=1/(1+exp(expVal));
end


pi0_W120=zeros(size(totU_W120));
piVal_W120=zeros(size(totU_W120 ));
for i=1:totU_W120
    [lPt,mPt,hPt]=fixed_point_v3(u_W120(i),theta);
    pi0_W120(i) = sum(xhist_W120{i}>mPt)/sum(xhist_W120{i}>=0);
    
    func=@(x) -g_func(x,u_W120(i),theta);
    delt_V_up=integral(func,hPt,mPt);
    delt_V_dwn=integral(func,lPt,mPt);

%     sig2_up=full(sigma2_func(hPt,u_W120(i),theta,120));
%     sig2_dwn=full(sigma2_func(lPt,u_W120(i),theta,120));
    
    sig2_up=full(sigma_EK_func(hPt,u_W120(i),theta,120));
    sig2_dwn=full(sigma_EK_func(lPt,u_W120(i),theta,120));

    V_dd_up=full(V_dd_func(hPt,u_W120(i),theta));
    V_dd_dwn=full(V_dd_func(lPt,u_W120(i),theta));

    expVal=2*(delt_V_dwn/sig2_dwn- delt_V_up/sig2_up)+(1/2)*(log(V_dd_up)-log(V_dd_dwn));
    %expVal=2*(delt_V_dwn/sig2_dwn - delt_V_up/sig2_up)+(1/2)*(log(V_dd_up_tst)-log(V_dd_dwn_tst));
    piVal_W120(i)=1/(1+exp(expVal));
end

c_est=glmfit(u_3,[round(pi0'*100) 100*ones(size(pi0'))],'binomial','logit');
c_est_W60=glmfit(u_W60,[round(pi0_W60'*100) 100*ones(size(pi0_W60'))],'binomial','logit');
c_est_W120=glmfit(u_W120,[round(pi0_W120'*100) 100*ones(size(pi0_W120'))],'binomial','logit');

c1_est=c_est(2)
c0_est=c_est(1)
cs_est=-c_est(1)/c_est(2);

c1_est_W60=c_est_W60(2)
c0_est_W60=c_est_W60(1)
cs_est_W60=-c_est_W60(1)/c_est_W60(2);

c1_est_W120=c_est_W120(2)
c0_est_W120=c_est_W120(1)
cs_est_W120=-c_est_W120(1)/c_est_W120(2);

logitFit = glmval(c_est,u_3,'logit');
logitFit_W60 = glmval(c_est_W60,u_W60,'logit');
logitFit_W120 = glmval(c_est_W120,u_W120,'logit');

subplot(2,1,1)
text(-.16,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
hold on
plot(u_W60(1:3:totU_W60),pi0_W60(1:3:totU_W60),'b+','LineWidth',2.0);
%plot(u_W60,piVal_W60,'b:','LineWidth',2.0);
plot(u_W60,logitFit_W60,'b')

plot(u_3,pi0,'k+','LineWidth',2.0);
%plot(u_3,piVal,'k:','LineWidth',2.0);
plot(u_3,logitFit,'k')

plot(u_W120,pi0_W120,'r+','LineWidth',2.0);
%plot(u_W120,piVal_W120,'r:','LineWidth',2.0);
plot(u_W120,logitFit_W120,'r')

xlim([0.1 .22])

hold off
xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
ylabel('\rho','FontSize',14,'FontWeight','bold','Color','k')
legend('Observed \rho, \Omega=60','Fit \rho, \Omega=60',...
    'Observed \rho, \Omega=90','Fit \rho, \Omega=90',...
        'Observed \rho, \Omega=120','Fit \rho, \Omega=120',...
            'Location','northwest')

saveas(f2,'Logist_Dip.png')

function b=bmtest(x)

    m3=skewness(x);
    m4=kurtosis(x);
    n=length(x);
    b=(m3^2+1)/(m4+3*( (n-1)^2/( (n-2)*(n-3) ) ) );

end

function val=g_func(x0,u,theta)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    val=a0+a*((u+x0).^n)./(K+(u+x0).^n)-x0;

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


function [dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)

    %  function		[dip,p_value,xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
    %
    % calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF (vector of sample values)
    % This routine calls the matlab routine 'HartigansDipTest' that actually calculates the DIP
    % NBOOT is the user-supplied sample size of boot-strap
    % Code by F. Mechler (27 August 2002)

    % calculate the DIP statistic from the empirical pdf
    [dip,xlow,xup, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf);
    N=length(xpdf);

    % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
    boot_dip=[];
    for i=1:nboot
       unifpdfboot=sort(unifrnd(0,1,1,N));
       [unif_dip]=HartigansDipTest(unifpdfboot);
       boot_dip=[boot_dip; unif_dip];
    end;
    boot_dip=sort(boot_dip);
    p_value=sum(dip<boot_dip)/nboot;

    % % Plot Boot-strap sample and the DIP statistic of the empirical pdf
    % figure(1); clf;
    % [hy,hx]=hist(boot_dip); 
    % bar(hx,hy,'k'); hold on;
    % plot([dip dip],[0 max(hy)*1.1],'r:');
end


function	[dip,xl,xu, ifault, gcm, lcm, mn, mj] = HartigansDipTest(xpdf)

    % function	[dip,xl,xu, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf)
    %
    % This is a direct translation by F. Mechler (August 27 2002)
    % into MATLAB from the original FORTRAN code of Hartigan's Subroutine DIPTST algorithm 
    % Ref: Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-325
    %
    % Appended by F. Mechler (September 2 2002) to deal with a perfectly unimodal input
    % This check the original Hartigan algorithm omitted, which leads to an infinite cycle
    %
    % HartigansDipTest, like DIPTST, does the dip calculation for an ordered vector XPDF using
    % the greatest convex minorant (gcm) and the least concave majorant (lcm),
    % skipping through the data using the change points of these distributions.
    % It returns the 'DIP' statistic, and 7 more optional results, which include
    % the modal interval (XL,XU), ann error flag IFAULT (>0 flags an error)
    % as well as the minorant and majorant fits GCM, LCM, and the corresponding support indices MN, and MJ

    % sort X in increasing order in column vector
    x=sort(xpdf(:));
    N=length(x);
    mn=zeros(size(x));
    mj=zeros(size(x));
    lcm=zeros(size(x));
    gcm=zeros(size(x));
    ifault=0;

    % Check that N is positive
    if (N<=0) 
       ifault=1;
       fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
       return;
    end;

    % Check if N is one
    if (N==1)
       xl=x(1);
       xu=x(N);
       dip=0.0;
       ifault=2;
       fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
       return;
    end;

    if (N>1)
       % Check that X is sorted
       if (x ~= sort(x))
          ifault=3;
          fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
          return;
       end;
       % Check for all values of X identical OR for case 1<N<4
       if ~((x(N)>x(1)) & (4<=N))
          xl=x(1);
          xu=x(N);
          dip=0.0;
          ifault=4;
          fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
          return;
       end;
    end;

    % Check if X is perfectly unimodal
    % Hartigan's original DIPTST algorithm did not check for this condition
    % and DIPTST runs into infinite cycle for a unimodal input
    % The condition that the input is unimodal is equivalent to having 
    % at most 1 sign change in the second derivative of the input p.d.f.
    xsign=-sign(diff(diff(x)));
    % This condition check below works even 
    % if the unimodal p.d.f. has its mode in the very first or last point of the input 
    % because then the boolean argument is Empty Matrix, and ANY returns 1 for an Empty Matrix
    posi=find(xsign>0);
    negi=find(xsign<0);
    if isempty(posi) | isempty(negi) | all(posi<min(negi))
       % A unimodal function is its own best unimodal approximation, with a zero corresponding dip
       xl=x(1);
       xu=x(N);
       dip=0.0;
       ifault=5;
        %fprintf(1,'\n  The input is a perfectly UNIMODAL input function\n');
       return;
    end;

    % LOW  contains the index of the current estimate of the lower end of the modal interval
    % HIGH contains the index of the current estimate of the upper end of the modal interval
    fn=N;
    low=1;
    high=N;
    dip=1/fn;
    xl=x(low);
    xu=x(high);

    % establish the indices over which combination is necessary for the convex minorant fit
    mn(1)=1;
    for j=2:N
       mn(j)=j-1;
       % here is the beginning of a while loop
       mnj=mn(j);
       mnmnj=mn(mnj);
       a=mnj-mnmnj;
       b=j-mnj;
       while ~( (mnj==1) | ((x(j)-x(mnj))*a < (x(mnj)-x(mnmnj))*b))
          mn(j)=mnmnj;
          mnj=mn(j);
          mnmnj=mn(mnj);
          a=mnj-mnmnj;
          b=j-mnj;
       end;   % here is the end of the while loop
    end; % end  for j=2:N

    % establish the indices over which combination is necessary for the concave majorant fit
    mj(N)=N;
    na=N-1;
    for jk=1:na
       k=N-jk;
       mj(k)=k+1;
       % here is the beginning of a while loop
       mjk=mj(k);
       mjmjk=mj(mjk);
       a=mjk-mjmjk;
       b=k-mjk;
       while ~( (mjk==N) | ((x(k)-x(mjk))*a < (x(mjk)-x(mjmjk))*b))
          mj(k)=mjmjk;
          mjk=mj(k);
          mjmjk=mj(mjk);
          a=mjk-mjmjk;
          b=k-mjk;
       end;   % here is the end of the while loop
    end; % end  for jk=1:na

    itarate_flag = 1;

    % start the cycling of great RECYCLE
    while itarate_flag 

    % collect the change points for the GCM from HIGH to LOW
    % CODE BREAK POINT 40
    ic=1;
    gcm(1)=high;
    igcm1=gcm(ic);
    ic=ic+1;
    gcm(ic)=mn(igcm1);
    while(gcm(ic) > low)
       igcm1=gcm(ic);
       ic=ic+1;
       gcm(ic)=mn(igcm1);
    end;
    icx=ic;

    % collect the change points for the LCM from LOW to HIGH
    ic=1;
    lcm(1)=low;
    lcm1=lcm(ic);
    ic=ic+1;
    lcm(ic)=mj(lcm1);
    while(lcm(ic) < high)
       lcm1=lcm(ic);
       ic=ic+1;
       lcm(ic)=mj(lcm1);
    end;
    icv=ic;

    % ICX, IX, IG are counters for the convex minorant
    % ICV, IV, IH are counters for the concave majorant
    ig=icx;
    ih=icv;

    % find the largest distance greater than 'DIP' between the GCM and the LCM from low to high
    ix=icx-1;
    iv=2;
    d=0.0;

    % Either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50;
    if ~(icx~=2 | icv~=2)
       d=1.0/fn;
    else
       iterate_BP50=1;
       while iterate_BP50
            % CODE BREAK POINT 50
            igcmx=gcm(ix);
          lcmiv=lcm(iv);
          if ~(igcmx > lcmiv)
             % if the next point of either the GCM or LCM is from the LCM then calculate distance here
             % OTHERWISE, GOTO BREAK POINT 55
             lcmiv1=lcm(iv-1);
             a=lcmiv-lcmiv1;
             b=igcmx-lcmiv1-1;
             dx=(x(igcmx)-x(lcmiv1))*a/(fn*(x(lcmiv)-x(lcmiv1)))-b/fn;
             ix=ix-1;
             if(dx < d) 
                goto60 = 1; 
             else
                d=dx;
                ig=ix+1;
                ih=iv;
                goto60 = 1;
             end;
          else
             % if the next point of either the GCM or LCM is from the GCM then calculate distance here
             % CODE BREAK POINT 55
             lcmiv=lcm(iv);
             igcm=gcm(ix);
             igcm1=gcm(ix+1);
             a=lcmiv-igcm1+1;
             b=igcm-igcm1;
             dx=a/fn-((x(lcmiv)-x(igcm1))*b)/(fn*(x(igcm)-x(igcm1)));
             iv=iv+1;
             if ~(dx < d) 
                d=dx;
                ig=ix+1;
                ih=iv-1;
             end;
             goto60 = 1;
          end;

          if goto60
             % CODE BREAK POINT 60
             if (ix < 1) ix=1; end;
             if (iv > icv) iv=icv; end;
             iterate_BP50 = (gcm(ix) ~= lcm(iv)); 
          end;
       end; % End of WHILE iterate_BP50
    end; % End of ELSE (IF ~(icx~=2 | icv~=2)) i.e., either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50

    % CODE BREAK POINT 65
    itarate_flag = ~(d < dip);
    if itarate_flag
    % if itarate_flag is true, then continue calculations and the great iteration cycle
    % if itarate_flag is NOT true, then stop calculations here, and break out of great iteration cycle to BREAK POINT 100

    % calculate the DIPs for the corrent LOW and HIGH

    % the DIP for the convex minorant
    dl=0.0;
    % if not true, go to CODE BREAK POINT 80
    if (ig ~= icx)
       icxa=icx-1;
       for j=ig:icxa
          temp=1.0/fn;
        jb=gcm(j+1);
          je=gcm(j);
          % if not true either, go to CODE BREAK POINT 74
          if ~(je-jb <= 1)
             if~(x(je)==x(jb))
                a=(je-jb);
                const=a/(fn*(x(je)-x(jb)));
                for jr=jb:je
                   b=jr-jb+1;
                   t=b/fn-(x(jr)-x(jb))*const;
                   if (t>temp) temp=t; end;
                end;
             end;
          end;
          % CODE BREAK POINT 74
          if (dl < temp) dl=temp; end;
       end;
    end;

    % the DIP for the concave majorant
    % CODE BREAK POINT 80
    du=0.0;
    % if not true, go to CODE BREAK POINT 90
    if ~(ih==icv)
       icva=icv-1;
       for k=ih:icva
          temp=1.0/fn;
          kb=lcm(k);
          ke=lcm(k+1);
          % if not true either, go to CODE BREAK POINT 86
          if ~(ke-kb <= 1)
             if ~(x(ke)==x(kb))
                a=ke-kb;
                const=a/(fn*(x(ke)-x(kb)));
                for kr=kb:ke
                   b=kr-kb-1;
                   t=(x(kr)-x(kb))*const-b/fn;
                   if (t>temp) temp=t; end;
                end;
             end;
          end;
          % CODE BREAK POINT 86
          if (du < temp) du=temp; end;
       end;
    end;

    % determine the current maximum
    % CODE BREAK POINT 90
    dipnew=dl;
    if (du > dl) dipnew=du; end;
    if (dip < dipnew) dip=dipnew; end;
    low=gcm(ig);
    high=lcm(ih);      

    end; % end of IF(itarate_flag) CODE from BREAK POINT 65

    % return to CODE BREAK POINT 40 or break out of great RECYCLE;
    end; % end of WHILE of great RECYCLE

    % CODE BREAK POINT 100
    dip=0.5*dip;
    xl=x(low);
    xu=x(high);

end


function u0=continAlg(theta,u0,dir,branch)
    ceq_val=[];
    import casadi.*
    x_sym = SX.sym('x_sym');
    u_sym = SX.sym('u_sym');

    a0_sym = SX.sym('a0_sym');
    a_sym = SX.sym('a_sym');
    K_sym = SX.sym('K_sym');
    n_sym = SX.sym('n_sym');

    theta_sym=[a0_sym a_sym K_sym n_sym];

    g=a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)-x_sym;
    %g_func = Function('g_func', {x_sym,u_sym,theta_sym}, {g});

    g_x=jacobian(g,x_sym);
    g_u=jacobian(g,u_sym);

    g_x_func = Function('g_x_func', {x_sym,u_sym,theta_sym}, {g_x});
    g_u_func = Function('g_u_func', {x_sym,u_sym,theta_sym}, {g_u});

    fold_loc=[0 0];
    flg=0;

    if branch==1
        [x0,~,~]=fixed_point_v3(u0,theta);
    else
        [~,~,x0]=fixed_point_v3(u0,theta);
    end
    stpSz=dir*.00001;
    while true

        gx=full(g_x_func(x0,u0,theta));
        gu=full(g_u_func(x0,u0,theta));
        A=[gx gu; 0 1];
        b=[0;1];
        if cond(A)>1e3
           break 
        end

        T=linsolve(A,b);

        x_hat=x0+stpSz*T(1);
        u_hat=u0+stpSz*T(2);

        mag=999;

        while abs(mag)>0.001
            gx=full(g_x_func(x_hat,u_hat,theta));
            gu=full(g_u_func(x_hat,u_hat,theta));
            gval=g_func(x_hat,u_hat,theta);

            A=-[gx gu; 0 1];
            b=[gval;0];
            if cond(A)>1e3
               flg=1;
               break 
            end

            T=linsolve(A,b);

            x_hat=x_hat+T(1);
            u_hat=u_hat+T(2);

            mag=max(abs(T));
        end

        if ~(abs(x0-x_hat)<0.5)||flg
            break
        else
            x0=x_hat;
            u0=u_hat;
        end
    end
    
end

