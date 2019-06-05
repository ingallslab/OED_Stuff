function FIM_vec=FIM_comp(pars,Omega,u_vec,bounds)

    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*

    theta=pars(1:4);
    c=pars(5:end);
    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);
    
    N_u=length(u_vec);
    Np=length(pars);
    N_c=length(c);
    N_th=length(theta);

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
    sigma2=-(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
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

    FIM_vec=zeros(Np,Np,N_u);

    for i=1:N_u

        if u_vec(i)<=bounds(1)
            %in lower monostable region
            [low,~,~]=fixed_point_v4(u_vec(i),theta);
            x_th=full(x_theta_func(low,u_vec(i),theta));
            sig2=full(sigma2_func(low,u_vec(i),theta,Omega));
            sig2_th=full(sigma2_theta_func(low,u_vec(i),theta,Omega));
            sig2_xval=full(sigma2_x_func(low,u_vec(i),theta,Omega));

            %FIM_vec(:,:,i)=(x_th'*x_th)./sig2+0.5*trace(sig2_th'*sig2_th)./(sig2^2);
            FIM_vec(1:N_th,1:N_th,i)=(x_th'*x_th)./sig2+0.5*((sig2_th+sig2_xval.*x_th)./sig2)'*((sig2_th+sig2_xval.*x_th)./sig2);

        elseif (bounds(1)<u_vec(i))&&(u_vec(i)<bounds(2))
            %in bistable region
            [low,~,high]=fixed_point_v4(u_vec(i),theta);
            std_low=sqrt(full(sigma2_func(low,u_vec(i),theta,Omega)));
            std_hgh=sqrt(full(sigma2_func(high,u_vec(i),theta,Omega)));
            
            lowBnd=min(low-4*std_low,high-4*std_hgh);
            upBnd=max(low+4*std_low,high+4*std_hgh);

            FIM_integrand=@(xo) (full(logLik_sens_func(xo,low,high,u_vec(i),pars,Omega))'*full(logLik_sens_func(xo,low,high,u_vec(i),pars,Omega)))*full(Lik_func(xo,low,high,u_vec(i),pars,Omega));
            %FIM_full = integral(FIM_integrand,-2e1,2e1,'ArrayValued',true);
            %should make bounds a function of the stnd dev,omega etc.
            FIM_full = integral(FIM_integrand,lowBnd,upBnd,'ArrayValued',true);

            FIM_vec(:,:,i)=FIM_full;

        else
            %in upper monostable region
            [~,~,high]=fixed_point_v4(u_vec(i),theta);
            x_th=full(x_theta_func(high,u_vec(i),theta));
            sig2=full(sigma2_func(high,u_vec(i),theta,Omega));
            sig2_th=full(sigma2_theta_func(high,u_vec(i),theta,Omega));
            sig2_xval=full(sigma2_x_func(high,u_vec(i),theta,Omega));

            %FIM_vec(:,:,i)=(x_th'*x_th)./sig2+0.5*trace(sig2_th'*sig2_th)./(sig2^2);
            FIM_vec(1:N_th,1:N_th,i)=(x_th'*x_th)./sig2+0.5*((sig2_th+sig2_xval.*x_th)./sig2)'*((sig2_th+sig2_xval.*x_th)./sig2);

        end
        FIM_vec(:,:,i)=FIM_vec(:,:,i).*(pars'*pars);
    end
