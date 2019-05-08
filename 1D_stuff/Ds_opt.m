
function [u_opt, w_opt]=Ds_opt(pars,Omega,u_vec,FIM_vec,bounds)
    
    addpath('/Users/nbraniff/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*

    N_u=length(u_vec);
    Np=length(pars);
    N_th=4;
    
    lam_vec={};
    lam0=[];
    Const={};
    lbw = [];lbg = [];
    ubw = [];ubg = [];
    lam_sum=0;
    FIM_tot=zeros(Np,Np);
    for i=1:N_u
        lambda_k = SX.sym(['lambda' '_' num2str(i)]);
        lbw = [lbw; 0]; ubw = [ubw,1];
        FIM_tot=FIM_tot+lambda_k*FIM_vec(:,:,i);

        lam_vec = {lam_vec{:}, lambda_k};
        lam0=[lam0; 1/N_u];
        lam_sum=lam_sum+lambda_k;
    end
    Const=[Const, {lam_sum-1}];
    lbg = [lbg; 0];ubg = [ubg; 0];

    %Is Ds optimal, rescaling invariant, (ie log par)
    J=-det(FIM_tot)./det(FIM_tot(N_th+1:end,N_th+1:end));

    % %% Create an NLP solver
    prob = struct('f', J, 'x', vertcat(lam_vec{:}), 'g', vertcat(Const{:}));
    % ipopt_opts=struct('linear_solver','ma27','hessian_approximation','exact',...
    %     'tol',1e-8,'dual_inf_tol',1,'constr_viol_tol',1e-4,...
    %     'acceptable_tol',1e-6,'acceptable_iter',15,'acceptable_dual_inf_tol',1e10,...
    %     'acceptable_obj_change_tol',1e20,'mu_strategy','adaptive');%,'max_iter',0);
    ipopt_opts=struct('linear_solver','ma27','hessian_approximation','exact',...
        'tol',1e-3,'dual_inf_tol',10,'constr_viol_tol',1e-2,...
        'acceptable_tol',1e-1,'acceptable_iter',8,'acceptable_dual_inf_tol',1e10,...
        'acceptable_obj_change_tol',1e10,'mu_strategy','adaptive','print_level',0,'max_iter',500);%,'max_iter',0);

    nlp_opts=struct('ipopt',ipopt_opts);
    solver = nlpsol('solver', 'ipopt', prob, nlp_opts);

    % Solve the NLP
    sol = solver('x0', lam0, 'lbx', lbw, 'ubx', ubw,...
                'lbg', lbg, 'ubg', ubg);
    lam = full(sol.x);
    
    lam(lam<1e-3)=0;
    boolLam=lam>0;
    nxtBool=logical(boolLam(1:end-1).*boolLam(2:end));

    flg=sum(nxtBool)~=0;

    cuttCnt=7;
    cntr=1;
    while flg&&cntr<cuttCnt
        new_u=[(u_vec([false; nxtBool])+u_vec([nxtBool; false]))/2];
        [u_vec,ind]=sort([u_vec new_u]);
        N_u=length(u_vec);
        FIM_tmp=cat(3,FIM_vec ,FIM_comp(pars,Omega,new_u,bounds));
        FIM_vec=FIM_tmp(:,:,ind);
        
        lam_vec={};
        lam0=[];
        Const={};
        lbw = [];lbg = [];
        ubw = [];ubg = [];
        lam_sum=0;
        FIM_tot=zeros(Np,Np);
        for i=1:N_u
            lambda_k = SX.sym(['lambda' '_' num2str(i)]);
            lbw = [lbw; 0]; ubw = [ubw,1];
            FIM_tot=FIM_tot+lambda_k*FIM_vec(:,:,i);

            lam_vec = {lam_vec{:}, lambda_k};
            lam0=[lam0; 1/N_u];
            lam_sum=lam_sum+lambda_k;
        end
        Const=[Const, {lam_sum-1}];
        lbg = [lbg; 0];ubg = [ubg; 0];

        %Is Ds optimal, rescaling invariant, (ie log par)
        J=-det(FIM_tot)./det(FIM_tot(N_th+1:end,N_th+1:end));

        % %% Create an NLP solver
        prob = struct('f', J, 'x', vertcat(lam_vec{:}), 'g', vertcat(Const{:}));
        nlp_opts=struct('ipopt',ipopt_opts);
        solver = nlpsol('solver', 'ipopt', prob, nlp_opts);

        % Solve the NLP
        sol = solver('x0', lam0, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg);
        lam = full(sol.x);
        
        lam(lam<1e-3)=0;
        boolLam=lam>0;
        nxtBool=logical(boolLam(1:end-1).*boolLam(2:end));
        
        flg=sum(nxtBool)~=0;
        cntr=cntr+1;
    end
    
    w_opt=lam(boolLam);
    u_opt=u_vec(boolLam);
    
    