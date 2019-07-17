addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

OMEGA=.01;

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


phi_1=[0.45 0];
phi_2=[0 45];

%%

x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);
g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;


%%
u_vec=linspace(0,1,200);
u_l=NaN; u_r=NaN;
flag=0;
sgn_cnt=1;
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

numRuns=10000;

num_bfr=20;
swtch_tm=200;

fn_tim=1200; %20000;
ss_tim=10000;
num_aft=round(fn_tim/30);

col=jet(num_bfr);

t_vec0=[0 logspace(0,log10(swtch_tm),num_bfr-1)];
%t_vec=[0 logspace(0,log10(fn_tim),num_aft-1)];
t_vec=[0 linspace(0,fn_tim,num_aft-1)];

g_roots_left=@(x_1) g_func1(x_1,phi_1(1),phi_1(2));
[x1_low,~,~]=fixed_point_v5(g_roots_left,numRuns);
x2_low=x2_null(x1_low,phi_1(1));

g_roots_right=@(x_1) g_func1(x_1,phi_2(1),phi_2(2));
[~,~,x1_high]=fixed_point_v5(g_roots_right,numRuns);
x2_high=x2_null(x1_high,phi_2(1));

%6000 W0.008 u.45 n500
%>10000 W0.02 u.45 n500
tmp=linspace(u_l,u_r,7);
Omeg_vec=[0.008 0.01 0.02];
%u_bist_vec=tmp(2:end-1);
u_bist_vec=[linspace(0.01,u_l-.05,6) linspace(u_l,u_r,12) linspace(u_r+.05,1,6)];
%u_bist_vec=[linspace(u_l,u_r,10)];

rws=4;
clms=6;

figs=zeros(size(Omeg_vec));

for l=1:length(Omeg_vec)
    Omeg_val=Omeg_vec(l);

    cur_dist_left=repmat([x1_low x2_low],numRuns,1);
    cur_dist_right=repmat([x1_high x2_high],numRuns,1);
    
    dist_left=zeros([size(cur_dist_left) num_bfr]);
    dist_right=zeros([size(cur_dist_right) num_bfr]);
    dist_left(:,:,1)=cur_dist_left;
    dist_right(:,:,1)=cur_dist_right;

    cnt=1;
    for i=2:num_bfr

        cur_dist_left0=Toggle_Ensemble(cur_dist_left ,phi_1,Omeg_val,t_vec0(i)-t_vec0(i-1),numRuns);
        cur_dist_right0=Toggle_Ensemble(cur_dist_right,phi_2,Omeg_val,t_vec0(i)-t_vec0(i-1),numRuns);

        dist_left(:,:,i)=cur_dist_left0;
        dist_right(:,:,i)=cur_dist_right0;   
    end

    stat1_left=zeros(1,num_bfr-1);
    stat2_left=zeros(1,num_bfr-1);
    stat1_right=zeros(1,num_bfr-1);
    stat2_right=zeros(1,num_bfr-1);
    for k=1:num_bfr-1
        [~,~,s1_left]= kstest2(dist_left(:,1,k),dist_left(:,1,end));
        [~,~,s2_left]= kstest2(dist_left(:,2,k),dist_left(:,2,end));

        [~,~,s1_right]= kstest2(dist_right(:,1,k),dist_right(:,1,end));
        [~,~,s2_right]= kstest2(dist_right(:,2,k),dist_right(:,2,end));

        stat1_left(k)=s1_left;
        stat2_left(k)=s2_left;

        stat1_right(k)=s1_right;
        stat2_right(k)=s2_right;
    end

    figs(l)=figure;
    slp_fits=zeros(2,length(u_bist_vec));
    leak_l=zeros(1,length(u_bist_vec));
    leak_r=zeros(1,length(u_bist_vec));
    logist_l=zeros(1,length(u_bist_vec));
    logist_r=zeros(1,length(u_bist_vec));
    dip_stat_l=zeros(1,length(u_bist_vec));
    dip_stat_r=zeros(1,length(u_bist_vec));
    %fig2=figure;
    for j=1:length(u_bist_vec)
        
        swtch_u=u_bist_vec(j);
        strng=['_W' num2str(Omeg_val) '_u' num2str(swtch_u)];

        u_vals=(phi_2-phi_1)*swtch_u+phi_1;

        cur_dist_left=cur_dist_left0;
        cur_dist_right=cur_dist_right0;
        
        dist_left=zeros([size(cur_dist_left) num_aft]);
        dist_right=zeros([size(cur_dist_right) num_aft]);
        dist_left(:,:,1)=cur_dist_left;
        dist_right(:,:,1)=cur_dist_right;

        for i=2:num_aft
            [l j i]

            cur_dist_left=Toggle_Ensemble(cur_dist_left ,u_vals,Omeg_val,t_vec(i)-t_vec(i-1),numRuns);
            cur_dist_right=Toggle_Ensemble(cur_dist_right,u_vals,Omeg_val,t_vec(i)-t_vec(i-1),numRuns);

            dist_left(:,:,i)=cur_dist_left;
            dist_right(:,:,i)=cur_dist_right;   
        end

        dist_left_all{l,j}=dist_left;
        dist_right_all{l,j}=dist_right;

        stat1_left=zeros(1,num_aft-1);
        stat2_left=zeros(1,num_aft-1);
        stat1_right=zeros(1,num_aft-1);
        stat2_right=zeros(1,num_aft-1);
        stat1_btwn=zeros(1,num_aft-1);
        stat2_btwn=zeros(1,num_aft-1);
              
        cnt=1;
        for k=1:num_aft-1

%             [~,~,s1_left]= kstest2(dist_left(:,1,k),dist_left(:,1,end));
%             [~,~,s2_left]= kstest2(dist_left(:,2,k),dist_left(:,2,end));
% 
%             [~,~,s1_right]= kstest2(dist_right(:,1,k),dist_right(:,1,end));
%             [~,~,s2_right]= kstest2(dist_right(:,2,k),dist_right(:,2,end));

            [~,~,s1_btwn]= kstest2(dist_left(:,1,k),dist_right(:,1,k));
            [~,~,s2_btwn]= kstest2(dist_left(:,2,k),dist_right(:,2,k));

%             stat1_left(k)=s1_left;
%             stat2_left(k)=s2_left;
% 
%             stat1_right(k)=s1_right;
%             stat2_right(k)=s2_right;

            stat1_btwn(k)=s1_btwn;
            stat2_btwn(k)=s2_btwn;
        end
       
%         subplot(2,2,1)
%         title('log-KS vs Time')
%         hold on
%         semilogy(t_vec(2:end),stat1_left,'b-');
%         semilogy(t_vec(2:end),stat2_left,'b.-');
%         semilogy(t_vec(2:end),stat1_right,'r-');
%         semilogy(t_vec(2:end),stat2_right,'r.-');
%         hold off
%         legend('X1-Left','X2-Left','X1-Right','X2-Right')
%         subplot(2,2,2)
%         title('log-KS vs Time')
%         hold on
%         semilogy(t_vec(2:end),stat1_btwn,'k-');
%         semilogy(t_vec(2:end),stat2_btwn,'k.-');
%         hold off
%         legend('X1-Btwn','X2-Btwn')
%         subplot(2,2,3)
%         title('log-KS vs Time')
%         hold on
%         semilogy(t_vec(2:end),stat1_left,'b-');
%         semilogy(t_vec(2:end),stat2_left,'b.-');
%         semilogy(t_vec(2:end),stat1_right,'r-');
%         semilogy(t_vec(2:end),stat2_right,'r.-');
%         hold off
%         xlim([0 1200])
%         ylim([.01 1])
%         legend('X1-Left','X2-Left','X1-Right','X2-Right')
        subplot(rws,clms,j)
        title('log-KS vs Time')
        hold on
        semilogy(t_vec(2:end),stat1_btwn,'k-');
        semilogy(t_vec(2:end),stat2_btwn,'k.-');
        hold off
        legend('X1-Btwn','X2-Btwn')
        xlim([0 1200])
        ylim([.1 1])
        
        t_sub=t_vec(2:end);
        inrng1=and(stat1_btwn<0.95,stat1_btwn>0.4);
        inrng2=and(stat2_btwn<0.95,stat2_btwn>0.4);
        p1 = polyfit(t_sub(inrng1),stat1_btwn(inrng1),1);
        p2 = polyfit(t_sub(inrng2),stat2_btwn(inrng2),1);
        slp_fits(1,j)=p1(1);
        slp_fits(2,j)=p2(1);
        
        if swtch_u<=u_l
            logist_l(j)=0;
            logist_r(j)=0;
            leak_l(j)=0;
            leak_r(j)=0;
        elseif swtch_u>=u_r
            logist_l(j)=1;
            logist_r(j)=1;
            leak_l(j)=1;
            leak_r(j)=1;
        else
            
            g_roots1=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
            [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
            x2_low=x2_null(x1_low,u_vals(1));
            x2_mid=x2_null(x1_mid,u_vals(1));
            x2_high=x2_null(x1_high,u_vals(1));
            
            J0=full(g_x_func([x1_mid x2_mid],u_vals(1),u_vals(2),theta));
            [eVecs,eVals]=eig(J0);
            eVals=diag(eVals)';
            vectr=eVecs(:,eVals<0);
            ortho_vectr=null(vectr');
            
            leak_l(j)=mean(((cur_dist_left-[x1_mid x2_mid])* -ortho_vectr)>0);
            leak_r(j)=mean(((cur_dist_right-[x1_mid x2_mid])* -ortho_vectr)>0);
            
            cur_dist_left=Toggle_Ensemble(cur_dist_left ,u_vals,Omeg_val,ss_tim-t_vec(end),numRuns);
            cur_dist_right=Toggle_Ensemble(cur_dist_right,u_vals,Omeg_val,ss_tim-t_vec(end),numRuns);
            
            coeff_left = pca(cur_dist_left);
            coeff_right = pca(cur_dist_right);
            
            dip_stat_l(j)=HartigansDipTest(cur_dist_left*coeff_left(:,1));
            dip_stat_r(j)=HartigansDipTest(cur_dist_right*coeff_left(:,1));
            
            logist_l(j)=mean(((cur_dist_left-[x1_mid x2_mid])* -ortho_vectr)>0);
            logist_r(j)=mean(((cur_dist_right-[x1_mid x2_mid])* -ortho_vectr)>0);

        end
        
        %savefig(fig4,['/Users/nbraniff/Documents/MATLAB/FOSBE_2019/2D_Figs/Hist_KS' strng '.fig']);
    end
    figure
    subplot(2,2,1)
    hold on
    plot(u_bist_vec,slp_fits(1,:))
    plot(u_bist_vec,slp_fits(2,:))
    hold off
    title('log-KS Slope Fits')
    
    subplot(2,2,2)
    hold on
    plot(u_bist_vec,leak_l)
    plot(u_bist_vec,leak_r)
    hold off
    title('Leaking at 20hrs')
    
    subplot(2,2,4)
    hold on
    plot(u_bist_vec,logist_l)
    plot(u_bist_vec,logist_r)
    hold off
    title('Eq. Logistic')
    
    subplot(2,2,3)
    hold on
    plot(u_bist_vec,dip_stat_l)
    plot(u_bist_vec,dip_stat_r)
    hold off
    title('Dip Statistic')
end

%save('/Users/nbraniff/Documents/MATLAB/FOSBE_2019/2D_Figs/dists','dist_left_all','dist_right_all','t_vec','Omeg_vec','u_bist_vec');

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
% figure
% hold on
% plot(u_tau,tau_low)
% plot(u_tau,tau_high)
% hold off

%%

% numRuns=5000;
% 
% strt_u=0;
% swtch_u=.48;
% 
% num_bfr=100;
% num_aft=100;
% num_tot=num_bfr+num_aft+1;
% swtch_tm=10;
% fn_tim=20000;
% 
% t_vec=[0 logspace(0,log10(swtch_tm),num_bfr), logspace(log10(swtch_tm+swtch_tm*.05),log10(fn_tim),num_aft)];
% u_vec=[strt_u repmat(strt_u,1,num_bfr), repmat(swtch_u,1,num_aft)];
% col=jet(num_tot);
% 
% u_vals=(phi_2-phi_1)*u_vec(1)+phi_1;
% g_roots=@(x_1) g_func1(x_1,u_vals(1),u_vals(2));
% [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots,numRuns);
% x2_low=x2_null(x1_low,u_vals(1));
% cur_dist=repmat([x1_low x2_low],numRuns,1);
% 
% fig1=figure;
% fig2=figure;
% dist=zeros([size(cur_dist) num_tot]);
% dist(:,:,1)=cur_dist;
% cnt=1;
% for i=2:num_tot
%     i
%     
%     u_vals=(phi_2-phi_1)*u_vec(i)+phi_1;
%     cur_dist=Toggle_Ensemble(cur_dist,u_vals,OMEGA,t_vec(i)-t_vec(i-1),numRuns);
%     
%     [f1,xi1] = ksdensity(cur_dist(:,1));
%     [f2,xi2] = ksdensity(cur_dist(:,2));
%     
%     dist(:,:,i)=cur_dist;
% 
%     if mod(i-1,round(num_tot/10))==0
%         figure(fig1)
%         subplot(1,2,1)
%         hold on
%         plot(xi1,f1,'Color',col(i,:));
%         hold off
%         subplot(1,2,2)
%         hold on
%         plot(xi2,f2,'Color',col(i,:));
%         hold off
%         
%         figure(fig2)
%         subplot(10,2,1+2*(cnt-1))
%         histogram(cur_dist(:,1),=-1:10:4000);
%         subplot(10,2,2+2*(cnt-1))
%         histogram(cur_dist(:,2),-1:10:3000);
%         cnt=cnt+1;
%     end
% 
% end
% 
% stat1=zeros(1,num_tot-1);
% stat2=zeros(1,num_tot-1);
% for k=1:num_tot-1
%     [~,~,s1]= kstest2(dist(:,1,k),dist(:,1,end));
%     [~,~,s2]= kstest2(dist(:,2,k),dist(:,2,end));
% 
%     stat1(k)=s1;
%     stat2(k)=s2;
% end
% figure
% plot(log10(t_vec(2:end)),stat1);
% figure
% plot(log10(t_vec(2:end)),stat2);


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
