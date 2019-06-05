addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%TO DO
%DOUBLE CHECK Ds OPTIMAL FORM, CAN WE USE IT WITH REGULAR D-opt??
%Log sensitivities
%How wide a quadrature range do we need, confirm with integral of Lik Func


%% Parameter declarations etc.
u_vec=[0:0.01:0.3];
N_u=length(u_vec);
bounds=[0.1 0.2];

%parameters
a0=0.5;
a=3;
K=9;
n=3;
Omega=90; %add this in

theta=[a0 a K n];
N_th=4;

bndryVal=8;
c0=bndryVal*(bounds(1)+bounds(2))/(bounds(1)-bounds(2));
c1=(-bndryVal-c0)/bounds(1);
c=[c0 c1];
N_c=2;

pars=[theta c];
Np=length(pars);

%%

%%%%cs_vec=[bounds(1)+.2*(diff(bounds)) bounds(1)+.5*(diff(bounds)) bounds(1)+.8*(diff(bounds))];
% cs_vec=[.1522 .1622 .1722];
% c1=111;
% N_c=length(cs_vec);
% 
% c0_vec=-cs_vec*c1;
% 
% clr=lines(2*N_c);
% 
% f1=figure;
% %f2=figure;
% for jj=1:N_c
%     
%     c0=-cs_vec(jj)*c1;
%     c=[c0 c1];
%     pars=[theta c];
%     
%     FIM_vec=FIM_comp(pars,Omega,u_vec,bounds);
%     
%     [u_Dopt, w_Dopt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds);
%     [u_Dsopt, w_Dsopt] = Ds_opt(pars,Omega,u_vec,FIM_vec,bounds);
%       
%     ys=1./(1+exp(-(c0+c1*u_vec)));
%     
%     figure(f1)
%     %set(gcf, 'InvertHardcopy', 'off')
%     %subplot(4,2,[1 3])
%     subplot(4,7,[1 2 3 8 9 10])
%     text(-.35,1.2,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     plot(u_vec,ys,'Color',clr(jj,:),'LineWidth',2.0);
%     hold off
%     %xlabel('Input, u','FontSize',18,'FontWeight','bold','Color','k')
%     ylabel('Branch Probability, \rho','FontSize',14,'FontWeight','bold','Color','k')
%     %title('c_0 Design Sensitivity','FontSize',20,'FontWeight','bold','Color','k')
%     if jj==N_c
%         legend(cellstr(num2str(c0_vec', 'c_0=%-.1f'))','Location','northwest')
%     end
%     
%     %subplot(4,2,5)
%     subplot(4,7,[15 16 17])
%     text(-.35,1.42,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     stem(u_Dsopt,w_Dsopt,'Color',clr(jj,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     hold off
%     %xlabel('Input, u','FontSize',18,'FontWeight','bold','Color','k')
%     ylabel({'D_s Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D_s Optimal Design','FontSize',20,'FontWeight','bold','Color','k')
%     
%     %subplot(4,2,7)
%     subplot(4,7,[22 23 24])
%     hold on
%     stem(u_Dopt,w_Dopt,'Color',clr(jj,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     text(-.35,1.42,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold off
%     xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel({'D Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D Optimal Design','FontSize',20,'FontWeight','bold','Color','k')
%     
% end
% 
% 
% c1_vec=[76,111,156];
% cs=.1622;
% N_c=length(c1_vec);
% 
% 
% for jj=1:N_c
%     
%     c1=c1_vec(jj);
%     c0=-cs*c1;
%     c=[c0 c1];
%     pars=[theta c];
%     
%     FIM_vec=FIM_comp(pars,Omega,u_vec,bounds);
%     
%     [u_Dopt, w_Dopt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds);
%     [u_Dsopt, w_Dsopt] = Ds_opt(pars,Omega,u_vec,FIM_vec,bounds);
%     
%     ys=1./(1+exp(-(c0+c1*u_vec)));
%     
%     %subplot(4,2,[2 4])
%     subplot(4,7,[5 6 7 12 13 14])
%     text(-.35,1.2,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     plot(u_vec,ys,'Color',clr(jj+3,:),'LineWidth',2.0);
%     hold off
%     %xlabel('Input, u','FontSize',18,'FontWeight','bold','Color','k')
%     ylabel('Branch Probability, \rho','FontSize',14,'FontWeight','bold','Color','k')
%     %title('c_1 Design Sensitivity','FontSize',20,'FontWeight','bold','Color','k')
%     if jj==N_c
%         legend(cellstr(num2str(c1_vec', 'c_1=%-.1f')),'Location','northwest')
%     end
%     
%     %subplot(4,2,6)
%     subplot(4,7,[19 20 21])
%     text(-.35,1.42,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     stem(u_Dsopt,w_Dsopt,'Color',clr(jj+3,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     hold off
%     %xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel({'D_s Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D_s Optimal Design','FontSize',20,'FontWeight','bold','Color','k')
% 
%     %subplot(4,2,8)
%     subplot(4,7,[26 27 28])
%     text(-.35,1.42,'F','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     stem(u_Dopt,w_Dopt,'Color',clr(jj+3,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     hold off
%     xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel({'D Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D Optimal Design','FontSize',20,'FontWeight','bold','Color','k')
% 
% end
% 
% saveas(f1,'LogistOED.png')

% runTot=10;
% runCnt=1;
% theta_vec=zeros(runTot,N_th);
% while runCnt<=runTot
%     
%     theta_drw=normrnd(theta,abs(.1*theta));
%     
%     a0=theta_drw(1);
%     a=theta_drw(2);
%     K=theta_drw(3);
%     n=theta_drw(4);
%     
%     tol=1e-4;
%     x_tst=[0:tol:6];
%     epsil=0.005;
%     g_bnd0_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0+x_tst).^theta_drw(4))./(theta_drw(3)+(0+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd1_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(1)-epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(1)-epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd2_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(2)+epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(2)+epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd3_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0.3+x_tst).^theta_drw(4))./(theta_drw(3)+(0.3+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     
%     if (g_bnd1_cnt==3)&&(g_bnd2_cnt==3)&&(g_bnd0_cnt==1)&&(g_bnd3_cnt==1)
%         theta_vec(runCnt,:)=theta_drw;
%         runCnt=runCnt+1;
%     end
%     
% end
% 
% cs=.1622;
% c1=111;
% c0=-cs*c1;
% c=[c0 c1];
% 
% htmap_D=zeros(N_c,N_u);
% htmap_Ds=zeros(N_c,N_u);
% 
% clr=lines(runTot);
% 
% f2=figure;
% for jj=1:runTot
% 
%     theta=theta_vec(jj,:);
%     pars=[theta c];
%     
%     FIM_vec=FIM_comp(pars,Omega,u_vec,bounds);
%     
%     [u_Dopt, w_Dopt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds);
%     [u_Dsopt, w_Dsopt] = Ds_opt(pars,Omega,u_vec,FIM_vec,bounds);
% 
%     u_mono1=[]; u_mono2=[];u_low=[];u_high=[];
%     roots_mono1=[]; roots_mono2=[];roots_low=[];roots_high=[];
%     flg=0;
%     for i=1:length(u_vec)
%         [lPt,mPt,hPt]=fixed_point_v3(u_vec(i),theta);
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
% 
%     figure(f2)
%     %set(gcf, 'InvertHardcopy', 'off')
%     subplot(4,1,1:2)
%     text(-.16,1.2,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     plot([u_mono1 u_low],[roots_mono1 roots_low],'Color',clr(jj,:),'LineWidth',2.0);
%     plot([u_high u_mono2],[roots_high roots_mono2],'Color',clr(jj,:),'LineWidth',2.0);
% %     plot(u_mono1,roots_mono1,'Color',clr(jj,:));
% %     plot(u_mono2,roots_mono2,'Color',clr(jj,:));
% %     plot(u_low,roots_low,'Color',clr(jj,:));
% %     plot(u_high,roots_high,'Color',clr(jj,:));
%     hold off
%     %xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel('Steady State, X','FontSize',14,'FontWeight','bold','Color','k')
%     %title('Bifurcation Plots for Sample of \theta','FontSize',20,'FontWeight','bold','Color','k')
%     
%     subplot(4,1,3)
%     text(-.16,1.4,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     stem(u_Dsopt,w_Dsopt,'Color',clr(jj,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     hold off
%     %xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel({'D_s Opt.';' Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D_s-opt Design','FontSize',20,'FontWeight','bold','Color','k')
% 
%     subplot(4,1,4)
%     text(-.16,1.4,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%     hold on
%     stem(u_Dopt,w_Dopt,'Color',clr(jj,:),'LineWidth',2.0)
%     xlim([0 0.3])
%     hold off
%     xlabel('Input, u','FontSize',14,'FontWeight','bold','Color','k')
%     ylabel({'D Opt.';' Weights'},'FontSize',14,'FontWeight','bold','Color','k')
%     %title('D-opt Design','FontSize',20,'FontWeight','bold','Color','k')
%     
% end
% 
% 
% saveas(f2,'ThetaOED.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% runTot=50;
% runCnt_D=1;
% runCnt_Ds=1;
% Ds_midpoints={};
% D_midpoints={};
% 
% cs=.1622;
% c1=111;
% c0=-cs*c1;
% c=[c0 c1];
% 
% f3=figure
% %set(gcf, 'InvertHardcopy', 'off')
% theta_vec=zeros(runTot,N_th);
% while runCnt_Ds<=runTot||runCnt_D<=runTot
%     
%     theta_drw=normrnd(theta,abs(.1*theta));
%     
%     a0=theta_drw(1);
%     a=theta_drw(2);
%     K=theta_drw(3);
%     n=theta_drw(4);
%     
%     tol=1e-4;
%     x_tst=[0:tol:6];
%     epsil=0.005;
%     g_bnd0_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0+x_tst).^theta_drw(4))./(theta_drw(3)+(0+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd1_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(1)-epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(1)-epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd2_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(2)+epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(2)+epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     g_bnd3_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0.3+x_tst).^theta_drw(4))./(theta_drw(3)+(0.3+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
%     
%     if (g_bnd1_cnt==3)&&(g_bnd2_cnt==3)&&(g_bnd0_cnt==1)&&(g_bnd3_cnt==1)
%         theta=theta_drw;
%         a=bounds(1)+.25*(diff(bounds)); 
%         b=bounds(1)+.75*(diff(bounds));
%         cs = a + (b-a).*rand;
%         c1 = 91 + (141-91).*rand;
%         c0=-cs*c1;
%         %c0=-(16 + (19-16).*rand);
%         c=[c0 c1];
%         pars=[theta c];
% 
%         FIM_vec=FIM_comp(pars,Omega,u_vec,bounds);
%         pars
%         if runCnt_D<=runTot
%             [u_Dopt, w_Dopt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds);
%             %if length(u_Dopt)==4
%                 D_midpoints{runCnt_D}=1./(1+exp(-(c0+c1* u_Dopt( and(bounds(1)<u_Dopt,u_Dopt<bounds(2)) ) )))';
%                 runCnt_D=runCnt_D+1
%             %end
%         end
%         if runCnt_Ds<=runTot
%             [u_Dsopt, w_Dsopt] = Ds_opt(pars,Omega,u_vec,FIM_vec,bounds);
%             %if length(u_Dsopt)==3
%                 Ds_midpoints{runCnt_Ds}=1./(1+exp(-(c0+c1*u_Dsopt(and(bounds(1)<u_Dsopt,u_Dsopt<bounds(2))))))';
%                 runCnt_Ds=runCnt_Ds+1
%             %end
%         end
%         if ((mod(runCnt_Ds,10)==0)&&(runCnt_Ds<=runTot))||((mod(runCnt_D,10)==0)&&(runCnt_D<=runTot))
%             subplot(2,1,1)
%             set(gca, 'FontSize', 20)
%             %histogram(Ds_midpoints,15,'Normalization','pdf')
%             histogram(vertcat(Ds_midpoints{:}),linspace(0,1,40))%,'Normalization','pdf')
%             text(-.16,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%             xlim([0 1]);
%             xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
%             ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
%             title('D_s Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')
% 
%             subplot(2,1,2)
%             set(gca, 'FontSize', 20)
%             %histogram(reshape(D_midpoints,2*length(D_midpoints),1),50,'Normalization','pdf')
%             histogram(vertcat(D_midpoints{:}),linspace(0,1,40))%,'Normalization','pdf')
%             text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
%             xlim([0 1]);
%             xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
%             ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
%             title('D Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')
% 
%             saveas(f3,'HistOED.png')
%         end
%     end
%     
% 
%     
% end
% 
% subplot(2,1,1)
% set(gca, 'FontSize', 20)
% %histogram(Ds_midpoints,15,'Normalization','pdf')
% histogram(vertcat(Ds_midpoints{:}),linspace(0,1,40))%,'Normalization','pdf')
% text(-.16,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
% xlim([0 1]);
% xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
% ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
% title('D_s Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')
% 
% subplot(2,1,2)
% set(gca, 'FontSize', 20)
% %histogram(reshape(D_midpoints,2*length(D_midpoints),1),50,'Normalization','pdf')
% histogram(vertcat(D_midpoints{:}),linspace(0,1,40))%,'Normalization','pdf')
% text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
% xlim([0 1]);
% xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
% ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
% title('D Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')
% 
% saveas(f3,'HistOED.png')

%% Hack

runTot=500;
runCnt_D=1;
runCnt_Ds=1;
Ds_midpoints={};
D_midpoints={};

cs=.1622;
c1=111;
c0=-cs*c1;
c=[c0 c1];

f4=figure
%set(gcf, 'InvertHardcopy', 'off')
theta_vec=zeros(runTot,N_th);
while runCnt_Ds<=runTot||runCnt_D<=runTot
    
    theta_drw=normrnd(theta,abs(.1*theta));
    
    a0=theta_drw(1);
    a=theta_drw(2);
    K=theta_drw(3);
    n=theta_drw(4);
    
    tol=1e-4;
    x_tst=[0:tol:6];
    epsil=0.005;
    g_bnd0_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0+x_tst).^theta_drw(4))./(theta_drw(3)+(0+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
    g_bnd1_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(1)-epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(1)-epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
    g_bnd2_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((bounds(2)+epsil+x_tst).^theta_drw(4))./(theta_drw(3)+(bounds(2)+epsil+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
    g_bnd3_cnt=sum(abs(diff(sign((theta_drw(1)+theta_drw(2).*((0.3+x_tst).^theta_drw(4))./(theta_drw(3)+(0.3+x_tst).^theta_drw(4))-x_tst)+eps)./2))>0);
    
    if (g_bnd1_cnt==3)&&(g_bnd2_cnt==3)&&(g_bnd0_cnt==1)&&(g_bnd3_cnt==1)
        theta=theta_drw;
        a=bounds(1)+.25*(diff(bounds)); 
        b=bounds(1)+.75*(diff(bounds));
        cs = a + (b-a).*rand;
        c1 = 91 + (141-91).*rand;
        c0=-cs*c1;
        %c0=-(16 + (19-16).*rand);
        c=[c0 c1];
        pars=[theta c];

        FIM_vec=FIM_comp(pars,Omega,u_vec,bounds);
        pars
        if runCnt_D<=runTot
            [u_Dopt, w_Dopt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds);
            %if length(u_Dopt)==4
                D_midpoints{runCnt_D}=1./(1+exp(-(c0+c1* u_Dopt( and(bounds(1)<u_Dopt,u_Dopt<bounds(2)) ) )))';
                runCnt_D=runCnt_D+1
            %end
        end
        if runCnt_Ds<=runTot
            [u_Dsopt, w_Dsopt] = Ds_opt(pars,Omega,u_vec,FIM_vec,bounds);
            %if length(u_Dsopt)==3
                Ds_midpoints{runCnt_Ds}=1./(1+exp(-(c0+c1*u_Dsopt(and(bounds(1)<u_Dsopt,u_Dsopt<bounds(2))))))';
                runCnt_Ds=runCnt_Ds+1
            %end
        end
        if ((mod(runCnt_Ds,10)==0)&&(runCnt_Ds<=runTot))||((mod(runCnt_D,10)==0)&&(runCnt_D<=runTot))
            subplot(2,1,1)
            set(gca, 'FontSize', 20)
            %histogram(Ds_midpoints,15,'Normalization','pdf')
            histogram(vertcat(Ds_midpoints{:}),linspace(0,1,100))%,'Normalization','pdf')
            text(-.16,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
            xlim([0 1]);
            xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
            ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
            title('D_s Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')

            subplot(2,1,2)
            set(gca, 'FontSize', 20)
            %histogram(reshape(D_midpoints,2*length(D_midpoints),1),50,'Normalization','pdf')
            histogram(vertcat(D_midpoints{:}),linspace(0,1,100))%,'Normalization','pdf')
            text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
            xlim([0 1]);
            xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
            ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
            title('D Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')

            saveas(f4,'HistOED2.png')
        end
    end
    

    
end

subplot(2,1,1)
set(gca, 'FontSize', 20)
%histogram(Ds_midpoints,15,'Normalization','pdf')
histogram(vertcat(Ds_midpoints{:}),linspace(0,1,100))%,'Normalization','pdf')
text(-.16,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
xlim([0 1]);
xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
title('D_s Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')

subplot(2,1,2)
set(gca, 'FontSize', 20)
%histogram(reshape(D_midpoints,2*length(D_midpoints),1),50,'Normalization','pdf')
histogram(vertcat(D_midpoints{:}),linspace(0,1,100))%,'Normalization','pdf')
text(-.16,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',26,'FontWeight','bold')
xlim([0 1]);
xlabel('\rho(u_{opt})','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Frequency','FontSize',14,'FontWeight','bold','Color','k')
title('D Optimal Sampling Percentiles','FontSize',18,'FontWeight','bold','Color','k')

saveas(f4,'HistOED2.png')




%% Functions

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

function val=g_func(x0,u,theta)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    val=a0+a*((u+x0).^n)./(K+(u+x0).^n)-x0;

end

