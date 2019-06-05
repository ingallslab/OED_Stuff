classdef OpExD_Model
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %deterministic model components
        theta       %model parameters
        u           %experimental inputs (input vars or time)
        mu          %symbolic expression for mu (scale and location parameters of the pdf)
        mu_func     %Casadi function for mu, call as mu_func(theta,u)
        
        %random sample/error density model
        x           %observation variable(s)
        pdf         %symbolic expression for the observation probability density function 
        pdf_func    %Casadi function for pdf, call as pdf_func(x,theta,u)
               
        
        %fitting and OED functions
        loglik      %symbolic expression for the model logliklihood
        loglik_func	%Casadi function logliklihood, call as loglik_func(x,theta,u)
        fim         %symbolix expression for model's expected fisher information matrix (FIM)
        fim_func    %Casadi function for FIM, called as fim(theta,u)
                
    end
    
    methods
        
        function obj = OpExD_Model(varargin)
          import casadi.*
          if nargin >= 4
             %set parameters (theta) and inputs (u)
             obj.theta=varargin{1};
             obj.u=varargin{2};
             
             %set mu symbolic and function
             obj.mu=varargin{3};
             obj.mu_func=Function('mu_func', {obj.theta,obj.u}, {obj.mu});
             switch varargin{4}
                 case 'Norm1_ConstVar'
                    %set obs vars (1D)
                    obj.x = SX.sym('x');
                    %set logliklihood symbolics and function
                    obj.pdf= 1/sqrt(2*pi*obj.mu(2)) * exp(-(obj.x-obj.mu(1))^2./(2*obj.mu(2)));
                    obj.pdf_func = Function('pdf_func', {obj.x,obj.theta,obj.u}, {obj.pdf});
                    %set logliklihood symbolics and function
                    obj.loglik = log(obj.pdf);%-0.5*log(2*pi*obj.mu(2)) - (obj.x-obj.mu(1))^2/(2*obj.mu(2));
                    obj.loglik_func = Function('loglik', {obj.x,obj.theta,obj.u}, {obj.loglik});
                    %set FIM symbolics and function
                    mu_theta=jacobian(obj.mu,obj.theta);
                    obj.fim=(mu_theta(1,:)'*mu_theta(1,:))./obj.mu(2);
                    obj.fim_func=Function('fim', {obj.theta,obj.u}, {obj.fim});
                 case 'Norm1_MeanVar'
                    %set obs vars (1D)
                    obj.x = SX.sym('x');
                    %set logliklihood symbolics and function
                    obj.pdf= 1/sqrt(2*pi*obj.mu(2)) * exp(-(obj.x-obj.mu(1))^2./(2*obj.mu(2)));
                    obj.pdf_func = Function('pdf_func', {obj.x,obj.theta,obj.u}, {obj.pdf});
                    %set logliklihood symbolics and function
                    obj.loglik = log(obj.pdf);%-0.5*log(2*pi*obj.mu(2)) - (obj.x-obj.mu(1))^2/(2*obj.mu(2));
                    obj.loglik_func = Function('loglik', {obj.x,obj.theta,obj.u}, {obj.loglik});
                    %set FIM symbolics and function
                    mu_theta=jacobian(obj.mu,obj.theta);
                    obj.fim=(mu_theta(1,:)'*mu_theta(1,:))./obj.mu(2)+(mu_theta(2,:)'*mu_theta(2,:))./obj.mu(2)^2;
                    obj.fim_func=Function('fim', {obj.theta,obj.u}, {obj.fim});
                case 'Norm1_UknwnVar'
                    %NEEDED
                case 'Log1_'
                    %NEEDED
                case 'Norm1_UknwnVar'
                    %NEEDED
                case 'Norm1_UknwnVar'
                    %NEEDED
                case 'Norm1_UknwnVar'
                    %NEEDED    
                case 'Norm1_UknwnVar'
                    %NEEDED
                case 'Norm2_UkwnVar'
                    %NEEDED
                case 'Norm2_KnwnVar'
                    %NEEDED
                case 'Custom'
                    obj.x = varargin{5};
                    p=varargin{6};
                    obj.pdf = Function('pdf', {obj.x,obj.theta,obj.u}, {p});
                otherwise
                    error('Unexpected distribution type. No PDF created.')
             end             
          else
              error('Improper number of inputs')
          end
        end
        
        
        function roe = optRelax(varargin)
          import casadi.*
          import approxOptExp.*
          if nargin >= 4
            obj=varargin{1};
            theta_in=varargin{2};
            u_grid=varargin{3};
            type=varargin{4};
             
            %settup problem for call to IPOPT, via Casadi
            Nu=length(u_grid);
            Np=length(theta_in);
            
            %create sampleing weight related vars
            xi_vec={};
            xi0=repmat(1/Nu,Nu,1);
            xi_sum=0;
            
            %constraint vars
            Const={};
            lbw = [];lbg = [];
            ubw = [];ubg = [];
            
            FIM_tot=zeros(Np,Np);
            FIM_init=zeros(Np,Np);
            for i=1:Nu
                xi_k = SX.sym(['xi' '_' num2str(i)]);
                lbw = [lbw; 0]; ubw = [ubw,1];
                xi_vec = {xi_vec{:}, xi_k};
                
                FIM_tot=FIM_tot+xi_k*obj.fim_func(theta_in,u_grid(i));
                FIM_init=FIM_init+(1/Nu)*full(obj.fim_func(theta_in,u_grid(i)));
                xi_sum=xi_sum+xi_k;
            end
            Const=[Const, {xi_sum-1}];
            lbg = [lbg; 0];ubg = [ubg; 0];

            switch type
                %The 'usual suspects', most common criteria
                case 'D'
                    J=-log(det(FIM_tot));
                case 'A'
                    J=trace(inv(FIM_tot));
                case 'E' 
                    %CVX PSD constraint method, from boyd convex prog
                    %book/examples
                    eig_min = SX.sym(['xi' '_' num2str(i)]);
                    lbw = [lbw; 0]; ubw = [ubw,Inf];
                    xi_vec = {xi_vec{:}, eig_min};
                    
                    %Using a constraint ot enforce PSD of FIM-eig*I, found
                    %in a paper https://vanderbei.princeton.edu/ps/sdp_nlp2.pdf
                    [D,L]=ldl(FIM_tot-eig_min*eye(Np,Np));
                    eig_min0=min(eig(FIM_init));
                    xi0=[xi0; eig_min0];
                    
                    for j=1:Np
                       Const=[Const, {D(j)}];
                       lbg = [lbg; 0];ubg = [ubg; Inf];
                    end
                    
                    J=-eig_min;
                    
                    %POWER ITERATION
                    %A=inv(FIM_tot);
                    %b=rand(Np,1);
                    %bs=b/norm(b);
                    
                    %for j=1:100
                    %    bs=A*bs;
                    %    bs=bs/norm(bs);
                    %end
                    
                    %J=norm(bs);
                    
                %criteria with additional inputs
                case 'Ds'
                    %Is Ds optimal, rescaling invariant, (ie log par)??????
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    N_targ=varargin{5};
                    J=-log(det(FIM_tot)./det(FIM_tot(N_targ+1:end,N_targ+1:end)));
                
                %prediction variance criteria
                case 'V'
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    u_targ=varargin{5};
                    Ntarg=length(u_targ);
                    
                    mu_theta=jacobian(obj.mu,obj.theta);
                    mu_theta_func=Function('mu_theta_func',{obj.theta,obj.u},{mu_theta});
                                        
                    L=zeros(Np,Np);
                    for i=1:Ntarg
                        L=L+mu_theta_func(theta_in,u_targ(i))'*mu_theta_func(theta_in,u_targ(i));
                    end
                    L=L/Ntarg;
                    
                    J=trace(L*inv(FIM_tot));
                
                case 'Custom'  
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    crit_func=varargin{5};
                    J=crit_func(FIM_tot);
                otherwise
                    error('Unexpected distribution type. No PDF created.')
            end

            % %% Create an NLP solver
            prob = struct('f', J, 'x', vertcat(xi_vec{:}), 'g', vertcat(Const{:}));
            ipopt_opts=struct('linear_solver','ma27','hessian_approximation','exact','mu_strategy','adaptive','max_iter',500,'print_level',0);

            nlp_opts=struct('ipopt',ipopt_opts,'print_time',0);
            solver = nlpsol('solver', 'ipopt', prob, nlp_opts);

            % Solve the NLP
            sol = solver('x0', xi0, 'lbx', lbw, 'ubx', ubw,...
                        'lbg', lbg, 'ubg', ubg);
                    
            status=solver.stats();
            display(status.return_status)
            xi = full(sol.x);
            
            if type=='E'
                xi=xi(1:end-1);
            end
            
            xi(xi<1e-3)=0;
            boolxi=xi>0;
            
            xi_opt=xi(boolxi)';
            u_opt=u_grid(boolxi);
            
            roe=rOptExp(u_opt,xi_opt);
             
          else
              error('Improper number of inputs')
          end
        end
        
        
        function eoe = optExact(varargin)
          import casadi.*
          import approxOptExp.*
          if nargin >= 4
            obj=varargin{1};
            theta_in=varargin{2};
            num_pts=varargin{3};
            u_strt=varargin{4};
            bnds=varargin{5};
            type=varargin{6};
             
            %settup problem for call to IPOPT, via Casadi
            Np=length(theta_in);
            
            u_vec={};
            lbu=bnds(1,:)';
            ubu=bnds(2,:)';
            dimU=length(obj.u);
            u0=reshape(u_strt,numel(u_strt),1);
           
            %constraint vars
            Const={};
            lbw = [];lbg = [];
            ubw = [];ubg = [];
            
            FIM_tot=zeros(Np,Np);
            for i=1:num_pts
                u_k = SX.sym(['xi' '_' num2str(i)],dimU);
                lbw = [lbw; lbu]; ubw = [ubw,ubu];
                u_vec = {u_vec{:}, u_k{:}};
                
                FIM_tot=FIM_tot+obj.fim_func(theta_in,u_k);
            end

            switch type
                %The 'usual suspects', most common criteria
                case 'D'
                    J=-log(det(FIM_tot));
                case 'A'
                    J=trace(inv(FIM_tot));
                case 'E' 
                    %CVX PSD constraint method, from boyd convex prog
                    %book/examples
                    eig_min = SX.sym(['xi' '_' num2str(i)]);
                    lbw = [lbw; 0]; ubw = [ubw,Inf];
                    xi_vec = {xi_vec{:}, eig_min};
                    
                    %Using a constraint ot enforce PSD of FIM-eig*I, found
                    %in a paper https://vanderbei.princeton.edu/ps/sdp_nlp2.pdf
                    [D,L]=ldl(FIM_tot-eig_min*eye(Np,Np));
                    eig_min0=min(eig(FIM_init));
                    xi0=[xi0; eig_min0];
                    
                    for j=1:Np
                       Const=[Const, {D(j)}];
                       lbg = [lbg; 0];ubg = [ubg; Inf];
                    end
                    
                    J=-eig_min;
                    
                    %POWER ITERATION
                    %A=inv(FIM_tot);
                    %b=rand(Np,1);
                    %bs=b/norm(b);
                    
                    %for j=1:100
                    %    bs=A*bs;
                    %    bs=bs/norm(bs);
                    %end
                    
                    %J=norm(bs);
                    
                %criteria with additional inputs
                case 'Ds'
                    %Is Ds optimal, rescaling invariant, (ie log par)??????
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    N_targ=varargin{5};
                    J=-log(det(FIM_tot)./det(FIM_tot(N_targ+1:end,N_targ+1:end)));
                
                %prediction variance criteria
                case 'V'
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    u_targ=varargin{5};
                    Ntarg=length(u_targ);
                    
                    mu_theta=jacobian(obj.mu,obj.theta);
                    mu_theta_func=Function('mu_theta_func',{obj.theta,obj.u},{mu_theta});
                                        
                    L=zeros(Np,Np);
                    for i=1:Ntarg
                        L=L+mu_theta_func(theta_in,u_targ(i))'*mu_theta_func(theta_in,u_targ(i));
                    end
                    L=L/Ntarg;
                    
                    J=trace(L*inv(FIM_tot));
                
                case 'Custom'  
                    if nargin~=5
                        error('Improper number of inputs')
                    end
                    crit_func=varargin{5};
                    J=crit_func(FIM_tot);
                otherwise
                    error('Unexpected distribution type. No PDF created.')
            end

            % %% Create an NLP solver
            prob = struct('f', J, 'x', vertcat(u_vec{:}), 'g', vertcat(Const{:}));
            ipopt_opts=struct('linear_solver','ma27','hessian_approximation','exact','mu_strategy','adaptive','max_iter',500,'print_level',0);

            nlp_opts=struct('ipopt',ipopt_opts,'print_time',0);
            solver = nlpsol('solver', 'ipopt', prob, nlp_opts);

            % Solve the NLP
            sol = solver('x0', u0, 'lbx', lbw, 'ubx', ubw,...
                        'lbg', lbg, 'ubg', ubg);
                    
            status=solver.stats();
            display(status.return_status)
            u_ext = full(sol.x);
                      
            eoe=eOptExp(u_ext);
             
          else
              error('Improper number of inputs')
          end
        end
        
        
        
        %Evaluat functions, for convenience
        function prob = evalMu(varargin)
          import casadi.*
          if nargin == 4
             obj=varargin{1};
             x_in=varargin{2};
             theta_in=varargin{3};
             u_in=varargin{4};
             prob=full(obj.mu_func(theta_in,u_in));
          else
              error('Improper number of inputs')
          end
        end
        
        function prob = evalPDF(varargin)
          import casadi.*
          if nargin == 4
             obj=varargin{1};
             x_in=varargin{2};
             theta_in=varargin{3};
             u_in=varargin{4};
             prob=full(obj.pdf_func(x_in,theta_in,u_in));
          else
              error('Improper number of inputs')
          end
        end
        
        function prob = evalLogLik(varargin)
          import casadi.*
          if nargin == 4
             obj=varargin{1};
             x_in=varargin{2};
             theta_in=varargin{3};
             u_in=varargin{4};
             prob=full(obj.loglik_func(x_in,theta_in,u_in));
          else
              error('Improper number of inputs')
          end
        end
        
        function prob = evalFIM(varargin)
          import casadi.*
          if nargin == 3
             obj=varargin{1};
             theta_in=varargin{2};
             u_in=varargin{3};
             prob=full(obj.fim_func(theta_in,u_in));
          else
              error('Improper number of inputs')
          end
        end

        
    end
    
end

