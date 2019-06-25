function syms = TwoD_Symbols_EM(uGrid1,uGrid2,numTrials)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all  
    
    u1_syms=cell(size(uGrid1,1),size(uGrid2,2));
    u2_syms=cell(size(uGrid1,1),size(uGrid2,2));
    x1_syms=cell(size(uGrid1,1),size(uGrid2,2),numTrials);
    x2_syms=cell(size(uGrid1,1),size(uGrid2,2),numTrials);
    for i=1:size(uGrid1,1)
        for j=1:size(uGrid2,2)
            u1_syms{i,j}=SX.sym(strcat('u1_syms_',num2str(i),'_',num2str(j)));
            u2_syms{i,j}=SX.sym(strcat('u2_syms_',num2str(i),'_',num2str(j)));
            for k=1:numTrials
                x1_syms{i,j,k} = SX.sym(strcat('x1_syms_',num2str(i),'_',num2str(j),'_',num2str(k)));
                x2_syms{i,j,k} = SX.sym(strcat('x2_syms_',num2str(i),'_',num2str(j),'_',num2str(k)));
            end
        end
    end
    
    c_11=SX.sym('c_11');
    c_12=SX.sym('c_12');
    c_21=SX.sym('c_21');
    c_22=SX.sym('c_22');
    
    detC = c_11*c_22-c_12*c_21;
    detC_func = Function('detC_func',{c_11,c_12,c_21,c_22},{detC});
    covar_sym=[c_11,c_12,c_21,c_22];
    
    x1_sym = SX.sym('x1_sym');
    x2_sym = SX.sym('x2_sym');
    x_sym=[x1_sym x2_sym];
    

    u1_sym = SX.sym('u1_sym');
    u2_sym = SX.sym('u2_sym');
    u_sym=[u1_sym u2_sym];
    
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

    
    g1_func = Function('g1_func',{theta_sym,u_sym,x_sym},{g1});
    g2_func = Function('g1_func',{theta_sym,u_sym,x_sym},{g2});
    
    g=[g1 g2];
    
    mu_1 = SX.sym('mu_1');
    mu_2 = SX.sym('mu_2');
    mu_sym = [mu_1 mu_2];
    
    invC_11 = c_22*detC^(-1);
    invC_22 = c_11*detC^(-1);
    invC_12 = -c_12*detC^(-1);
    invC_21 = -c_21*detC^(-1);
    invC = [invC_11 invC_12 invC_21 invC_22];
    
    qForm = invC_11*(mu_1 - x1_sym)^2 + invC_22*(mu_2-x2_sym)^2 + (invC_12+invC_21)*(mu_1-x1_sym)*(mu_2-x2_sym);
    qForm_func = Function('qForm',{mu_sym, covar_sym, x_sym},{qForm});
    gaussian = ((2*pi*detC)^(-1))*exp(-0.5*qForm);
    gauss_func = Function('gaussian',{mu_sym,covar_sym,x_sym},{gaussian});
    
    optimSyms = {alpha1_sym, alpha2_sym, beta1_sym, beta2_sym, K1_sym, K2_sym, n1_sym, n2_sym, kappa1_sym, kappa2_sym, m1_sym, m2_sym};
    
    lUB=[inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    lLB=[-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf];
    
    nlC={};
    nlUB=[];
    nlLB=[];
    
    QFunc=0;
    for i=1:size(u1_syms,1)
        for j=1:size(u2_syms,2)           
            xstar_top_1_ij = SX.sym(strcat('xstar_top_1_',num2str(i),'_',num2str(j)));
            optimSyms={optimSyms{:}, xstar_top_1_ij};
            lUB = [lUB; inf];
            lLB = [lLB; -inf];
            xstar_top_2_ij = SX.sym(strcat('xstar_top_2_',num2str(i),'_',num2str(j)));
            optimSyms={optimSyms{:}, xstar_top_2_ij};
            lUB = [lUB; inf];
            lLB = [lLB; -inf];
            xstar_bot_1_ij = SX.sym(strcat('xstar_bot_1_',num2str(i),'_',num2str(j)));
            optimSyms={optimSyms{:}, xstar_bot_1_ij};
            lUB = [lUB; inf];
            lLB = [lLB; -inf];
            xstar_bot_2_ij = SX.sym(strcat('xstar_bot_2_',num2str(i),'_',num2str(j)));
            optimSyms={optimSyms{:}, xstar_bot_2_ij};
            lUB = [lUB; inf];
            lLB = [lLB; -inf];
            
            covar_ij_top = [SX.sym(strcat('covarT_',num2str(i),'_',num2str(j),'_11'))...
                            SX.sym(strcat('covarT_',num2str(i),'_',num2str(j),'_12'))...
                            SX.sym(strcat('covarT_',num2str(i),'_',num2str(j),'_21'))...
                            SX.sym(strcat('covarT_',num2str(i),'_',num2str(j),'_22'))];
            logdetCij_top = detC_func(covar_ij_top(1),covar_ij_top(2),covar_ij_top(3),covar_ij_top(4));
            optimSyms={optimSyms{:}, covar_ij_top(1),covar_ij_top(2),covar_ij_top(3),covar_ij_top(4)};
            lUB = [lUB; inf; inf; inf; inf];
            lLB = [lLB; -inf; -inf; -inf; -inf];
            
            covar_ij_bot = [SX.sym(strcat('covarB_',num2str(i),'_',num2str(j),'_11'))...
                            SX.sym(strcat('covarB_',num2str(i),'_',num2str(j),'_12'))...
                            SX.sym(strcat('covarB_',num2str(i),'_',num2str(j),'_21'))...
                            SX.sym(strcat('covarB_',num2str(i),'_',num2str(j),'_22'))];
            logdetCij_bot = detC_func(covar_ij_bot(1),covar_ij_bot(2),covar_ij_bot(3),covar_ij_bot(4));
            optimSyms={optimSyms{:}, covar_ij_bot(1),covar_ij_bot(2),covar_ij_bot(3),covar_ij_bot(4)};
            lUB = [lUB; inf; inf; inf; inf];
            lLB = [lLB; -inf; -inf; -inf; -inf];
            
            w_ij_top=SX.sym(strcat('w_',num2str(i),'_',num2str(j),'_top'));
            optimSyms={optimSyms{:}, w_ij_top};
            lUB=[lUB;inf];
            lLB=[lLB;0];
            w_ij_bot=SX.sym(strcat('w_',num2str(i),'_',num2str(j),'_bot'));
            optimSyms={optimSyms{:}, w_ij_bot};
            lUB=[lUB;inf];
            lLB=[lLB;0];
            
            wSum = w_ij_top + w_ij_bot;
            nlC=[nlC, {wSum}];
            nlUB=[nlUB;1];
            nlLB=[nlLB;1];     
            
            gStar1 = g1_func(theta_sym, [u1_syms{i,j} u2_syms{i,j}], [xstar_top_1_ij xstar_bot_1_ij]);
            gStar2 = g2_func(theta_sym, [u1_syms{i,j} u2_syms{i,j}], [xstar_top_2_ij xstar_bot_2_ij]);
            
            nlC = [nlC, {gStar1}];
            nlUB = [nlUB; 0];
            nlLB = [nlLB; 0];
            nlC = [nlC, {gStar2}];
            nlUB = [nlUB; 0];
            nlLB = [nlLB; 0];
            
            for k=1:numTrials
                gau_top = gauss_func([xstar_top_1_ij xstar_top_2_ij], covar_ij_top, [x1_syms{i,j,k} x2_syms{i,j,k}]);
                gau_bot = gauss_func([xstar_bot_1_ij xstar_bot_2_ij], covar_ij_bot, [x1_syms{i,j,k} x2_syms{i,j,k}]);
                
                weight_top = w_ij_top*gau_top/(w_ij_top*gau_top+w_ij_bot*gau_bot);
                weight_bot = w_ij_bot*gau_bot/(w_ij_top*gau_top+w_ij_bot*gau_bot);
                
                TermOne = weight_top*(log(w_ij_top)-0.5*logdetCij_top)+weight_bot*(log(w_ij_bot)-0.5*logdetCij_bot);
                
                TermTwo = -0.5*...
                    (weight_top*(covar_ij_top(1)*(xstar_top_1_ij - x1_syms{i,j,k})^2 + covar_ij_top(4)*(xstar_top_2_ij-x2_syms{i,j,k})^2 +...
                    (covar_ij_top(2)+covar_ij_top(3))*(xstar_top_1_ij-x1_syms{i,j,k})*(xstar_top_2_ij-x2_syms{i,j,k}))+...
                    weight_bot*(covar_ij_bot(1)*(xstar_bot_1_ij - x1_syms{i,j,k})^2 + covar_ij_bot(4)*(xstar_bot_2_ij-x2_syms{i,j,k})^2 +...
                    (covar_ij_bot(2)+covar_ij_bot(3))*(xstar_bot_1_ij-x1_syms{i,j,k})*(xstar_bot_2_ij-x2_syms{i,j,k})));
                
                QFunc=QFunc+TermOne+TermTwo;
            end
        end
    end
    p = [];
    p = [p; vertcat(x1_syms{:})];
    p = [p; vertcat(x2_syms{:})];
    p = [p; vertcat(u1_syms{:})];
    p = [p; vertcat(u2_syms{:})];
    QFunc=-QFunc;
    q=Function('q',{vertcat(optimSyms{:}),p},{QFunc});
    nlp = struct('x', vertcat(optimSyms{:}), 'f', QFunc, 'g', vertcat(nlC{:}), 'p', p);
    syms = struct('nlp',nlp,'lbx',lLB,'ubx',lUB,'lbg',nlLB,'ubg',nlUB,'qfunc',q);
end