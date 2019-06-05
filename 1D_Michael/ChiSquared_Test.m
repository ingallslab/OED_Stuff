function g = ChiSquared_Test(theta)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*

    x_sym = SX.sym('x_sym');
    u_sym = SX.sym('u_sym');
    Omega_sym = SX.sym('Omega_sym');

    a0_sym = SX.sym('a0_sym');
    a_sym = SX.sym('a_sym');
    K_sym = SX.sym('K_sym');
    n_sym = SX.sym('n_sym');

    theta_sym=[a0_sym a_sym K_sym n_sym];
    theta_true = [0.5,3,9,3];

    c0_sym = SX.sym('c0_sym');
    c1_sym = SX.sym('c1_sym');

    c_sym=[c0_sym c1_sym];

    par_sym=[theta_sym c_sym];
    par_true = [0.5,3,9,3,-18.1,111.7];

    g=a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)-x_sym;
    %g_func = Function('g_func', {x_sym,u_sym,theta_sym}, {g});
    g_x=jacobian(g,x_sym);
    g_x_func = Function('g_x_func', {x_sym,u_sym,theta_sym}, {g_x});


    %sigma2=(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)^2/Omega_sym;
    sigma2=-(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
    %log sensitivities
    sigma2_theta=jacobian(sigma2,theta_sym);
    sigma2_x=jacobian(sigma2,x_sym);
    %sigma2_theta=jacobian(sigma2,theta_sym).*theta_sym;

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


    u_monostable1=[]; u_monostable2=[];u_lowbranch=[];u_highbranch=[];
    x_monostable1=[]; x_monostable2=[];x_lowbranch=[];x_highbranch=[];
    flg=0;
    uVals=linspace(0,0.3,100);
    for i=1:length(uVals)
        [lPt,~,hPt]=fixed_point_v4(uVals(i),theta_true);

        if (lPt==hPt&&flg==0)
            u_monostable1=[u_monostable1 uVals(i)];
            x_monostable1=[x_monostable1 lPt];

        elseif(lPt==hPt&&flg==1)
            u_monostable2=[u_monostable2 uVals(i)];
            x_monostable2=[x_monostable2 lPt];

        else
            u_lowbranch=[u_lowbranch uVals(i)];
            x_lowbranch=[x_lowbranch lPt];
            u_highbranch=[u_highbranch uVals(i)];
            x_highbranch=[x_highbranch hPt];

            flg=1;
        end 
    end


    
    g=0; 
    for i=1:length(u_monostable1)
        g=g+...
            integral(@(x) full(2*x*(log(phi_func(x,x_monostable1(i),u_monostable1(i),theta,90))-...
            log(phi_func(x,x_monostable1(i),u_monostable1(i),theta_true,90)))),0,4);
    end
    for i=1:length(u_monostable2)
        g=g+...
            integral(@(x) full(2*x*(log(phi_func(x,x_monostable2(i),u_monostable2(i),theta,90))-...
            log(phi_func(x,x_monostable2(i),u_monostable2(i),theta_true,90)))),0,4);
    end
    for i=1:length(u_highbranch)
        g=g+...
            integral(@(x) full(2*x*(logLik_func(x,x_lowbranch(i),x_highbranch(i),u_lowbranch(i),par_sym,90)-...
            logLik_func(x,x_lowbranch(i),x_highbranch(i),u_lowbranch(i),par_true,90))),0,4);
    end
end


