function likelihood = computeLikelihood_v3(xVals, uVals, params, sigma)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    
    likelihood = 0;
    
    a0_sym = SX.sym('a0');
    a_sym = SX.sym('a');
    K_sym = SX.sym('K');
    n_sym = SX.sym('n');
    theta_sym = [a0_sym, a_sym, K_sym, n_sym]
    
    x_sym = SX.sym('x');
    u_sym = SX.sym('u');
    sigma_sym = SX.sym('sigma');
    
    g_con = a0_sym + (a_sym*(x_sym+u_sym).^n_sym)./(K+(u_sym+x_sym).^n_sym) - x_sym;
    g_con_func = Function('g_con_func',{x_sym, u_sym, theta_sym}, {g_con});
    
    c0_sym = SX.sym('c0');
    c1_sym = SX.sym('c1');
    rho = 1/(1+exp(-(c0_sym + c1_sym*u_sym)));
    
    x_obs_sym = SX.sym('x_obs');
    gaussian = (1/(sigma_sym * sqrt(2*pi))) * exp(-((x_sym - x_obs_sym).^2)/(2*sigma_sym.^2));
    
    
end