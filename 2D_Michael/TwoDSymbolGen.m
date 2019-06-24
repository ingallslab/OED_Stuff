function syms = TwoDSymbolGen()
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all

    x1_sym = SX.sym('x1_sym');
    x2_sym = SX.sym('x2_sym');
    x_sym=[x1_sym x2_sym];

    u1_sym = SX.sym('u1_sym');
    u2_sym = SX.sym('u2_sym');
    Omega_sym = SX.sym('Omega_sym');

    alpha1_sym = SX.sym('alpha1_sym');
    alpha2_sym = SX.sym('alpha2_sym');
    beta1_sym = SX.sym('beta1_sym');
    beta2_sym = SX.sym('beta2_sym');
    K1_sym = SX.sym('K1_sym');
    K2_sym = SX.sym('K2_sym');
    n1_sym  = SX.sym('n1_sym');
    n2_sym = SX.sym('n2_sym');
    kappa1_sym = SX.sym('kappa1_sym');
    kappa2_sym = SX.sym('kappa2_sym');
    m1_sym = SX.sym('m1_sym');
    m2_sym = SX.sym('m2_sym');

    theta1_sym=[alpha1_sym beta1_sym K1_sym n1_sym kappa1_sym m1_sym];
    theta2_sym=[alpha2_sym beta2_sym K2_sym n2_sym kappa2_sym m2_sym];
    theta_sym=[alpha1_sym alpha2_sym beta1_sym beta2_sym K1_sym K2_sym n1_sym n2_sym kappa1_sym kappa2_sym m1_sym m2_sym];

    g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u2_sym./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
    g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u1_sym./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;

    g=[g1 g2];
    g_x=jacobian(g,x_sym);
    g_x_func = Function('g_x_func', {x_sym,u1_sym,u2_sym,theta_sym}, {g_x});

    %g_u_func = Function('g_u_func', {x_sym,u_sym,theta_sym}, {g_u});
    %g_theta=jacobian(g,theta_sym);
    %g1_x=jacobian(g1,x_sym);
    %g2_x=jacobian(g2,x_sym);
    %g_u=jacobian(g,u_sym);
    %full(g_x_func([500 300],u0_1,u0_2,theta))
    
end