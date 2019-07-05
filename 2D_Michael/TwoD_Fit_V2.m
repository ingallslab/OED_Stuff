%Brute force non Casadi approach to EM.

umax_1=1;
umax_2=100;
[u1_grid,u2_grid] = meshgrid(0:0.1:umax_1,0:10:umax_2);
numTrials=20;
SSAData=cell(size(u1_grid,1),size(u1_grid,2));
x1vals=[];
x2vals=[];
u1vals=[];
u2vals=[];
I=size(u1_grid,1);
J=size(u2_grid,2);
for i=1:size(u1_grid,1)
    for j=1:size(u2_grid,2)
        s = dlmread(strcat('Data/SSAData_',num2str(i),'_',num2str(j)),'\t');
        SSAData{i,j}=s;        
        x1vals=[x1vals; s(:,1)];
        x2vals=[x2vals; s(:,2)];
        u1vals=[u1vals;u1_grid(i,j)];
        u2vals=[u2vals;u2_grid(i,j)];
    end
end

covar_11=ones(I*J,1)'*100;
covar_21=zeros(I*J,1)'*100;
covar_22=ones(I*J,1)'*100;

w1=1/2;
w2=1/2;

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

theta_true=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];
x0={theta_true,covar_11,covar_21,covar_22,w1,w2};
x0=vertcat([x0{:}]);

obj = @(x)EW(x,u1vals,u2vals,x1vals,x2vals);

options=optimset('Display','iter');
FIT=fminsearch(obj,x0,options)
dlmwrite('Fit_Results.txt',FIT,'\t');

function E=EW(x,uVals1,uVals2,xVals1,xVals2)
    I=11;
    J=11;
    E=expectationStep(uVals1,uVals2,xVals1,xVals2,x(1:12),x(13),x(14),x(15:14+I*J),x(15+I*J:14+2*I*J),x(15+2*I*J:14+3*I*J),...
        x(15:14+I*J),x(15+I*J:14+2*I*J),x(15+2*I*J:14+3*I*J));
end

function ex=expectationStep(u1_grid,u2_grid,x1_vals,x2_vals,theta_guess,w1,w2,covar_1_11,covar_1_12,covar_1_22,covar_2_11,covar_2_12,covar_2_22)
    %u0_1=0.25;
    %u0_2=20;

    I=11;
    J=11;
    K=20;
    u1_grid=reshape(u1_grid,[I,J]);
    u2_grid=reshape(u2_grid,[I,J]);
    x1_vals=reshape(x1_vals,[K,J*I]);
    x2_vals=reshape(x2_vals,[K,J*I]);
    x1_vals=reshape(x1_vals,[K,J,I]);
    x2_vals=reshape(x2_vals,[K,J,I]);
    covar_1_11=reshape(covar_1_11,[I,J]);
    covar_1_12=reshape(covar_1_12,[I,J]);
    covar_1_22=reshape(covar_1_22,[I,J]);
    covar_2_11=reshape(covar_2_11,[I,J]);
    covar_2_12=reshape(covar_2_12,[I,J]);
    covar_2_22=reshape(covar_2_22,[I,J]);

    xyz_stable_1=[];
    xyz_stable_2=[];
    xyz_unstable_1=[];
    xyz_unstable_2=[];

    alpha_1 =theta_guess(1);
    alpha_2 =theta_guess(2);
    beta_1  =theta_guess(3);
    beta_2  =theta_guess(4);
    K_1     =theta_guess(5); 
    K_2     =theta_guess(6);
    n_1     =theta_guess(7);
    n_2     =theta_guess(8); 
    kappa_1 =theta_guess(9);
    kappa_2 =theta_guess(10);
    m_1     =theta_guess(11);
    m_2     =theta_guess(12);
    
    %x1_null=@(x_2,u_2) alpha_1 + beta_1./(1+((x_2./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1);
    x2_null=@(x_1,u_1) alpha_2 + beta_2./(1+((x_1./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2);

    g_func1=@(x_1,u_1,u_2) alpha_1 + beta_1*1./(1+(((alpha_2 +beta_2*1./(1+((x_1./K_1)*(1/(1+(u_1./kappa_1).^m_1))).^n_2))./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x_1;
    %g_func2=@(x_2,u_1,u_2) alpha_2 + beta_2*1./(1+(((alpha_1 +beta_1*1./(1+((x_2./K_2)*(1/(1+(u_2./kappa_2).^m_2))).^n_1))./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x_2;

    %g_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)-x(1) alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)-x(2)];

    %B_func=@(x,u_1,u_2) [alpha_1 + beta_1./(1+((x(2)./K_2)*(1./(1+(u_2./kappa_2).^m_2))).^n_1)+x(1) 0;...
    %    0 alpha_2 + beta_2./(1+((x(1)./K_1)*(1./(1+(u_1./kappa_1).^m_1))).^n_2)+x(2)];
   
    tmp=0;
    for i=1:size(u1_grid,1)
        %disp(strcat('i = ',num2str(i)));
        for j=1:size(u2_grid,2)
            co_1=[[covar_1_11(i,j) covar_1_12(i,j)];[covar_1_12(i,j) covar_1_22(i,j)]];
            co_2=[[covar_2_11(i,j) covar_2_12(i,j)];[covar_2_12(i,j) covar_2_22(i,j)]];
            
            g_roots1=@(x_1) g_func1(x_1,u1_grid(i,j),u2_grid(i,j));
            [x1_low,x1_mid,x1_high]=fixed_point_v5(g_roots1,3000);
            
            if (x1_low==x1_mid&&x1_mid==x1_high)
                xyz_stable_1=[xyz_stable_1; [u1_grid(i,j) u2_grid(i,j) x1_low]];
                
                x2_low=x2_null(x1_low,u1_grid(i,j));
                xyz_stable_2=[xyz_stable_2; [u1_grid(i,j) u2_grid(i,j) x2_low]];
               
                tmp=tmp-0.5*log10(det(co_1));
                for k=1:20
                    X=[x1_vals(k,i,j), x2_vals(k,i,j)]-[x1_low, x2_low];
                    tmp=tmp-0.5*X*(co_1\X');
                end
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
                tmp2=0;
                for k=1:20
                    %disp([[x1_vals(k,i,j), x2_vals(k,i,j)],[x1_high, x2_high]]);
                    weight_1=w1*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_high, x2_high],co_1)/...
                        (w1*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_high, x2_high],co_1)+...
                        w2*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_low, x2_low],co_2));
                    weight_2=w2*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_low, x2_low],co_2)/...
                        (w1*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_high, x2_high],co_1)+...
                        w2*mvnpdf([x1_vals(k,i,j), x2_vals(k,i,j)],[x1_low, x2_low],co_2));

                    tmp2=tmp2+weight_1*log10(w1)+weight_2*log10(w2);
                    tmp2=tmp2-0.5*weight_1*log10(det(co_1))-0.5*weight_2*log10(det(co_2));

                    X_1=[x1_vals(k,i,j), x2_vals(k,i,j)]-[x1_high, x2_high];
                    tmp2=tmp2-0.5*weight_1*X_1*(co_1\X_1');
                    X_2=[x1_vals(k,i,j), x2_vals(k,i,j)]-[x1_low, x2_low];
                    tmp2=tmp2-0.5*weight_2*X_2*(co_2\X_2');
                    tmp=tmp+tmp2;
                end
            end
            
        end
    end
    ex=-tmp;
end
