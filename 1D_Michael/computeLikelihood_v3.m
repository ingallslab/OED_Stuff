function likelihood = computeLikelihood_v3(xVals, uVals, params)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    
    Omega = 90;
    uLow = 0.1;
    uHigh = 0.2;
    
    theta = [params(1),params(2),params(3),params(4)];
    a0= theta(1);
    a = theta(2);
    K = theta(3);
    n = theta(4);
    Bounds = [0.1,0.2];
    
    tol=1e-4;
    g_bnd1_cnt=0;
    maxU=5;
    while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    g_bnd2_cnt=0;
    maxU=5;
    while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
        likelihood=1e8;
        return
    end
    
    
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
    g_x=jacobian(g,x_sym);
    
    sigma2=(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
    pi0=1/(1+exp(-(c0_sym+c1_sym*u_sym)));
    
    x_obs_sym = SX.sym('x_obs_sym');
    phi=(1/(sqrt(2*pi*sigma2)))*exp((-(x_obs_sym-x_sym).^2)./(2*sigma2));
    phi_func = Function('phi_func', {x_obs_sym,x_sym,u_sym,theta_sym,Omega_sym}, {phi});
    
    x_low_sym = SX.sym('x_low_sym');
    x_high_sym = SX.sym('x_high_sym');
    
    Lik=(1-pi0)*phi_func(x_obs_sym,x_low_sym,u_sym,theta_sym,Omega_sym)+pi0*phi_func(x_obs_sym,x_high_sym,u_sym,theta_sym,Omega_sym);
    
    logLik=log(Lik);
    logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_sym,par_sym,Omega_sym}, {logLik});
    
    logLik_tot = 0;
    
    xstars = []; %mean x values
    constraints=[]; %nonlinear constraint functions
    lbg = []; %nonlinear constraint bounds lower
    ubg = []; %nonlinear constraint bounds upper
    lbw = []; %linear lower bounds
    ubw = []; %linear upper bounds
    
    for i =1:length(uVals)
        u = uVals(i);
        
        constraints = [constraints {g}];
        lbg = [lbg; 0; 0];
        ubg = [ubg; 0; 0];
        for j = 1:size(xVals,1)
            if u<uLow
                xstar_i = SX.sym(strcat('xstar_',num2str(i)));
                xstars=[xstars xstars_i];
                logLik_tot = logLik_tot + phi_func(xVals(j,i),xstar_i,u,par_sym,Omega_sym);
            elseif u>uHigh
                xstar_i = SX.sym(strcat('xstar_',num2str(i)));
                xstars=[xstars xstars_i];
                logLik_tot = logLik_tot + phi_func(xVals(j,i),xstar_i,u,par_sym,Omega_sym);
            else
                
            end
        end
    end
    
    logLikelihood = Function('logLikelihood', {xstars,par_sym,Omega_sym},{logLik_tot});
    loglikelihood
    
    
end

