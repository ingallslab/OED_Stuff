function sym = generateOptimSymbols(xVals,uVals)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    Omega = 90;
    uLow = 0.1;
    uHigh = 0.2;   
    
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
    
    sigma2=-(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
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
    
    xvals = cell(size(xVals,1),size(xVals,2));
    for i=1:size(xVals,1)
        for j=1:size(xVals,2)
            x_symij = SX.sym(strcat('xvals_',num2str(i),num2str(j)));
            xvals{i,j}=x_symij;
        end
    end
    
    
    xstars_m = []; %mean x values -> monostable region
    constraints=[]; %nonlinear constraint functions
    lbg = []; %nonlinear constraint bounds lower
    ubg = []; %nonlinear constraint bounds upper
    lbw = []; %linear lower bounds
    ubw = []; %linear upper bounds
    
    xstars_h =[]; %mean x values -> upper branch, bistable region
    xstars_l =[]; %mean x values -> lower branch, bistable region
    disp('setting up symbolic expression');
    for i =1:length(uVals)
        u = uVals(i);
        
        if u<uLow
            xstar_i = SX.sym(strcat('xstar_',num2str(i)));
            xstars_m=[xstars_m; xstar_i];

            gstar=a0_sym+a_sym*((u+xstar_i).^n_sym)./(K_sym+(u+xstar_i).^n_sym)-xstar_i;
            constraints = [constraints, {gstar}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:size(xVals,1)
                logLik_tot = logLik_tot + log(phi_func(xvals{j,i},xstar_i,u,theta_sym,Omega));
            end
        elseif u>uHigh
            xstar_i = SX.sym(strcat('xstar_',num2str(i)));
            xstars_m=[xstars_m; xstar_i];

            gstar=a0_sym+a_sym*((u+xstar_i).^n_sym)./(K_sym+(u+xstar_i).^n_sym)-xstar_i;
            constraints = [constraints, {gstar}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:size(xVals,1)
                logLik_tot = logLik_tot + log(phi_func(xvals{j,i},xstar_i,u,theta_sym,Omega));
            end
        else
            xstar_h = SX.sym(strcat('xstar_h_',num2str(i)));
            xstars_h=[xstars_h; xstar_h];
            xstar_l = SX.sym(strcat('xstar_l_',num2str(i)));
            xstars_l=[xstars_l; xstar_l];

            constraints = [constraints, {xstar_l - xstar_h}];
            lbg = [lbg; -inf];
            ubg = [ubg; -0.001];

            gstar_h=a0_sym+a_sym*((u+xstar_h).^n_sym)./(K_sym+(u+xstar_h).^n_sym)-xstar_h;
            constraints = [constraints, {gstar_h}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            gstar_l=a0_sym+a_sym*((u+xstar_l).^n_sym)./(K_sym+(u+xstar_l).^n_sym)-xstar_l;
            constraints = [constraints, {gstar_l}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:size(xVals,1)
                logLik_tot = logLik_tot + logLik_func(xvals{j,i},xstar_l,xstar_h,u,par_sym,Omega);
            end
        end
    end
    cons = vertcat(constraints{:});
    p = vertcat(xvals{:});
    optimVars = [xstars_m; xstars_h; xstars_l; transpose(par_sym)];
    nlp = struct('x', optimVars, 'f', -logLik_tot, 'g', cons, 'p', p);
    options.error_on_fail = true;
    options.ipopt.max_iter = 40;
    solver = nlpsol('solver','ipopt',nlp,options);
    disp('solver has generated, beginning optimization');
    sym = struct('optimVars',optimVars, 'logLik_tot',-logLik_tot,'lbg',lbg,'ubg',ubg,'solver',solver,'fixedparams',p);
end