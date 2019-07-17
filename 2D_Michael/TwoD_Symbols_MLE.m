function syms = TwoD_Symbols_MLE(xVals,uVals)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    
    phi_1=[0.45 0];
    phi_2=[0 45];

    Phi=[phi_1;phi_2];
    
    x1_sym = SX.sym('x1_sym');
    x2_sym = SX.sym('x2_sym');
    x_sym=[x1_sym x2_sym];
    Omega_sym = SX.sym('Omega_sym');
    u_cntrl_sym = SX.sym('u_cntrl_sym');
    u_val_sym=(phi_2-phi_1)*u_cntrl_sym+phi_1;

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
    theta_sym=[alpha1_sym alpha2_sym beta1_sym beta2_sym K1_sym K2_sym n1_sym n2_sym kappa1_sym kappa2_sym m1_sym m2_sym];

    gamma0_sym = SX.sym('gamma0_sym');
    gamma1_sym = SX.sym('gamma1_sym');
    gamma_sym=[gamma0_sym gamma1_sym];

    par_sym=[theta_sym gamma_sym];

    x1_low_sym = SX.sym('x1_low_sym');
    x2_low_sym = SX.sym('x2_low_sym');
    x_low_sym=[x1_low_sym x2_low_sym];

    x1_high_sym = SX.sym('x1_high_sym');
    x2_high_sym = SX.sym('x2_high_sym');
    x_high_sym=[x1_high_sym x2_high_sym];

    x1_obs_sym = SX.sym('x1_obs_sym');
    x2_obs_sym = SX.sym('x2_obs_sym');
    x_obs_sym=[x1_obs_sym x2_obs_sym];

    g1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)-x1_sym;
    g2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)-x2_sym;
    g=[g1 g2];

    b1 = alpha1_sym + beta1_sym./(1+((x2_sym./K2_sym)*(1./(1+(u_val_sym(2)./kappa2_sym).^m2_sym))).^n1_sym)+x1_sym;
    b2 = alpha2_sym + beta2_sym./(1+((x1_sym./K1_sym)*(1./(1+(u_val_sym(1)./kappa1_sym).^m1_sym))).^n2_sym)+x2_sym;
    B = [b1 0; 0 b2];
    A=jacobian(g,x_sym);
    A_func = Function('A_func',{x_sym,u_cntrl_sym,theta_sym},{A});
    
    B_func = Function('B_func',{x_sym,u_cntrl_sym,theta_sym},{B});
    
    %?????????????????????????????????????????????????????????
    C_vec=(inv(kron(eye(size(A)),A)+kron(A,eye(size(A))))*vec(-B))./Omega_sym;
    Cvec_func = Function('Cvec_func',{x_sym,u_cntrl_sym,theta_sym,Omega_sym},{C_vec});
    
    
    %?????????????????????????????????????????????????????????
    C=reshape(C_vec,size(A)); 
    Cinv=inv(C);
    Cinv_func = Function('Cinv_func', {x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {Cinv});
    C_theta=jacobian(C_vec([1 2 4]),theta_sym);

    pdf=exp(-.5*(x_obs_sym-x_sym)*inv(C)*(x_obs_sym-x_sym)')/sqrt(det(C)*(2*pi)^2);
    pdf_func = Function('pdf_func', {x_obs_sym,x_sym,u_cntrl_sym,theta_sym,Omega_sym}, {pdf});

    pdf_low=pdf_func(x_obs_sym,[x1_low_sym,x2_high_sym],u_cntrl_sym,theta_sym,Omega_sym);
    pdf_high=pdf_func(x_obs_sym,[x1_high_sym,x2_low_sym],u_cntrl_sym,theta_sym,Omega_sym);

    rho=1/(1+exp(-(gamma1_sym*(u_cntrl_sym-gamma0_sym))));
    Lik=rho*pdf_low+(1-rho)*pdf_high;
    Lik_func = Function('Lik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {Lik});

    logLik=log(Lik);
    logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym}, {logLik});
    
    logLik_xH = gradient(logLik,x_high_sym);
    logLik_xH_func = Function('logLik_xH_func',{x_obs_sym,x_low_sym,x_high_sym,u_cntrl_sym,par_sym,Omega_sym},{logLik_xH});
    
    xvals = cell(size(xVals));
    for i=1:size(xvals)
        for j=1:size(xVals{i},1)
            xvals{i}=[xvals{i};[SX.sym(strcat('x1_',num2str(i),'_',num2str(j))) SX.sym(strcat('x2_',num2str(i),'_',num2str(j)))]];
        end
    end
    
    x1stars_h =[]; 
    x1stars_l =[];
    x1stars_m =[]; 
    x2stars_h =[]; 
    x2stars_l =[];
    x2stars_m =[]; 
    constraints=[];
    lbg=[]; %nonlinear constraint bounds (associated to $constraints)
    ubg=[];
    lbw=[]; %linear constraint bounds (associated to $optimVars)
    ubw=[];
    
    logLik_tot=0;
    
    disp('setting up symbolic expression');
    
    uLow = 18.0019;
    uHigh = 20.2516;
    for i =1:length(uVals)
        uT = uVals(i);
        uV=(phi_2-phi_1)*uT+phi_1;
        u1=uV(1);
        u2=uV(2);
        uvec=[uV(1),uV(2)];
        u=norm(uvec);
        if u<uLow
            x1star_i = SX.sym(strcat('x1star_',num2str(i)));
            x1stars_m=[x1stars_m; x1star_i];
            x2star_i = SX.sym(strcat('x2star_',num2str(i)));
            x2stars_m=[x2stars_m; x2star_i];
            
            g1s = alpha1_sym + beta1_sym./(1+((x2star_i./K2_sym)*(1./(1+(u2./kappa2_sym).^m2_sym))).^n1_sym)-x1star_i;
            g2s = alpha2_sym + beta2_sym./(1+((x1star_i./K1_sym)*(1./(1+(u1./kappa1_sym).^m1_sym))).^n2_sym)-x2star_i;
            
            lbw=[lbw;0];
            ubw=[ubw;inf];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            
            constraints = [constraints, {g1s}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            constraints = [constraints, {g2s}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            for j = 1:size(xvals{i},1)
                logLik_tot = logLik_tot + log(pdf_func(xvals{i}(j,:),[x1star_i x2star_i],u1,theta_sym,Omega_sym));
            end
        elseif u>uHigh
            x1star_i = SX.sym(strcat('x1star_',num2str(i)));
            x1stars_m=[x1stars_m; x1star_i];
            x2star_i = SX.sym(strcat('x2star_',num2str(i)));
            x2stars_m=[x2stars_m; x2star_i];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            g1s = alpha1_sym + beta1_sym./(1+((x2star_i./K2_sym)*(1./(1+(u2./kappa2_sym).^m2_sym))).^n1_sym)-x1star_i;
            g2s = alpha2_sym + beta2_sym./(1+((x1star_i./K1_sym)*(1./(1+(u1./kappa1_sym).^m1_sym))).^n2_sym)-x2star_i;

            constraints = [constraints, {g1s}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            constraints = [constraints, {g2s}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            for j = 1:size(xvals{i},1)
                logLik_tot = logLik_tot + log(pdf_func(xvals{i}(j,:),[x1star_i x2star_i],u1,theta_sym,Omega_sym));
            end
        else
            x1star_h = SX.sym(strcat('x1star_h_',num2str(i)));
            x1stars_h=[x1stars_h; x1star_h];
            x1star_l = SX.sym(strcat('x1star_l_',num2str(i)));
            x1stars_l=[x1stars_l; x1star_l];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            x2star_h = SX.sym(strcat('x2star_h_',num2str(i)));
            x2stars_h=[x2stars_h; x2star_h];
            x2star_l = SX.sym(strcat('x2star_l_',num2str(i)));
            x2stars_l=[x2stars_l; x2star_l];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            lbw=[lbw;0];
            ubw=[ubw;inf];
            
            g1sh = alpha1_sym + beta1_sym./(1+((x2star_h./K2_sym)*(1./(1+(u2./kappa2_sym).^m2_sym))).^n1_sym)-x1star_h;
            g2sh = alpha2_sym + beta2_sym./(1+((x1star_h./K1_sym)*(1./(1+(u1./kappa1_sym).^m1_sym))).^n2_sym)-x2star_h;
            
            g1sl = alpha1_sym + beta1_sym./(1+((x2star_l./K2_sym)*(1./(1+(u2./kappa2_sym).^m2_sym))).^n1_sym)-x1star_l;
            g2sl = alpha2_sym + beta2_sym./(1+((x1star_l./K1_sym)*(1./(1+(u1./kappa1_sym).^m1_sym))).^n2_sym)-x2star_l;
            
            constraints = [constraints, {x1star_l - x1star_h}];
            lbg = [lbg; -inf];
            ubg = [ubg; -0.001];
            constraints = [constraints, {x2star_l - x2star_h}];
            lbg = [lbg; 0.001];
            ubg = [ubg; inf];
            
            constraints = [constraints, {g1sh}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            constraints = [constraints, {g1sl}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            constraints = [constraints, {g2sh}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            constraints = [constraints, {g2sl}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
            for j = 1:size(xvals{i},1)
                logLik_tot = logLik_tot + logLik_func(xvals{i}(j,:),[x1star_l x2star_l],[x1star_h x2star_h],u1,par_sym,Omega_sym);
            end
        end
    end
    logLik_tot = -logLik_tot;
    
    lbw=[lbw;0.00001*ones(14,1)];
    ubw=[ubw;inf*ones(14,1)];
    
    disp('Symbols Done, Generating NLP');
    cons = vertcat(constraints{:});
    p = vertcat(xvals{:});
    p = vertcat(p(:));
    p = [p;Omega_sym];
    optimVars = [x1stars_m; x1stars_h; x1stars_l; x2stars_m; x2stars_h; x2stars_l; (par_sym)'];
    nlp = struct('x', optimVars, 'f', logLik_tot, 'g', cons, 'p', p);
    
    options.error_on_fail = true;
    options.ipopt.max_iter = 20;
    options.ipopt.mu_strategy = 'adaptive'; %Adaptive mu
    %options.monitor = {'nlp_f','nlp_g'}; %Allows monitoring of values
    solver = nlpsol('solver','ipopt',nlp,options);
    
    disp('solver has generated, beginning optimization');
    LF = Function('logLik_tot',{optimVars,p},{logLik_tot});
    lF = @(a,b) LF(a,b); %Numerical output
    
    %This just spits out the symbols in a struct
    syms = struct('theta_sym',theta_sym,'optimVars',optimVars, 'logLik_tot',logLik_tot,'loglikF',...
        lF,'cons',Function('cons',{optimVars,p},{cons}),...
        'lbg',lbg,'ubg',ubg,'solver',solver,'fixedparams',p,'lbw',lbw,'ubw',ubw,...
        'loglik_xH',logLik_xH_func,'lik_func',Lik_func,'cinv',Cinv_func,'cvec',Cvec_func,'Afunc',A_func,'Bfunc',B_func);
end
