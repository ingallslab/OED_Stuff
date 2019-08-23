classdef OpExD_Model
    %Primary class for OpExD package
    
    properties
        
        %deterministic model components
        theta       %model parameters
        u           %experimental inputs (input vars or time)
        mu          %symbolic expression for mu (scale and location parameters of the pdf)
        mu_func     %Casadi function for mu, call as mu_func{i}(theta,u)
        mu_vec_func %concatenates all the means from each pdf

        %random sample/error density model
        y           %observation variable(s)
        y_vec       %observation variables concatenated
        pdf         %symbolic expression for the observation probability density function of y, depends on mu, u
        pdf_func    %Casadi function for pdf, call as pdf_func{i}(x,theta,u)
        pdf_vec_func%concatenates the pdfs into a vector function
        
        %fitting and OED functions
        loglik      %symbolic expression for the model logliklihood
        loglik_func	%Casadi function logliklihood, call as loglik_func(x,theta,u)
        loglik_tot  %sum of all loglikelihoods from each pdf
        loglik_tot_func
        fim         %symbolix expression for model's expected fisher information matrix (FIM)
        fim_func    %Casadi function for FIM, called as fim(theta,u)
        fim_tot     %adds up all the FIM matrices from each pdf
        fim_tot_func
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Class Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = OpExD_Model(varargin)
            import casadi.*
            if nargin >= 1
                %models contains structs representing each model type
                % Example for two models: 
                % models={ struct('Type',{'Normal'}, 'Mean',{mu1_sym}, ...
                %           'Covariance',{cov1_sym}, 'Details',{'2x2 normpdf'},...
                %           'Params',{theta1_sym}, 'Input',{u1_sym}),
                %          struct('Type',{'Normal'}, 'Mean',{mu2_sym}, ...
                %           'Covariance',{cov2_sym}, 'Details',{'4z4 normpdf'},...
                %           'Params',{theta2_sym}, 'Input',{u2_sym})}
                %Here, sym_vec is a vector containing the individual symbols
                %required for the models - this allows us to build the loglik_func 
                %and fim_func objects. Details is a string describing the model for debugging purposes.
                %Params are the desired parameters we need to fit to.
                
                obj.mu_func = {};
                obj.pdf={};
                
                models    = varargin{1};
                obj.theta = varargin{2};
                obj.u     = varargin{3};
                
                obj.fim={};
                obj.fim_func={};
                
                obj.loglik={};
                obj.loglik_func={};
                
                for i=1:length(models)
                    dist=models{i}.Type;
                    switch dist
                        case 'Normal'
                            dim = length(models{i}.Mean);

                            obj.mu{i}=models{i}.Mean;
                            obj.mu_func{i}=Function('mu_func', {obj.theta,obj.u}, {models{i}.Mean});
                            
                            sigma2 = det(models{i}.Covariance);
                            invCo = inv(models{i}.Covariance);
                            
                            obj.y{i}=SX.sym('y',dim)';
                            
                            %set pdf symbolics and function
                            obj.pdf{i}= 1/sqrt(((2*pi)^dim)*sigma2) * ...
                                exp(-0.5*(obj.y{i}'-models{i}.Mean)'*invCo*(obj.y{i}'-models{i}.Mean));
                            obj.pdf_func{i} = Function('pdf_func', {obj.y{i},obj.theta,obj.u}, {obj.pdf{i}});
                            
                            %set logliklihood symbolics and function
                            obj.loglik{i} = log(obj.pdf{i}); 
                            obj.loglik_func{i} = Function('loglik', {obj.y{i},obj.theta,obj.u}, {obj.loglik{i}});
                            
                            %set FIM symbolics and function
                            mu_theta=jacobian(obj.mu{i},obj.theta);
                            k = length(obj.theta);
                            I = {};
                            for m=1:k
                                for n=1:k
                                    sig_m = reshape(jacobian(models{i}.Covariance,obj.theta(m)),size(models{i}.Covariance));
                                    sig_n = reshape(jacobian(models{i}.Covariance,obj.theta(n)),size(models{i}.Covariance));
                                    I{m,n} = mu_theta(:,m)'*invCo*mu_theta(:,n)+...
                                        0.5*trace(invCo*sig_m*invCo*sig_n);
                                end
                            end
                            I = reshape(vertcat(I{:}),k,k);
                            obj.fim{i}=I;
                            obj.fim_func{i}=Function('fim', {obj.theta,obj.u}, {obj.fim{i}});    
                        case 'Lognorm'
                            
                        case 'Pois'
                            
                        case 'CustomAll'
                            %user provides expressions for pdf as well as FIM (i.e.
                            %they take the expectation analytically
                            
                            
                        case 'CustomPDF'
                            %user provides expressions for pdf only, we numerically
                            %integrate it to get the FIM
                            x=models{i}.x;
                            p=models{i}.PDF;
                            obj.theta{i}=models{i}.theta;
                            obj.u{i}=models{i}.u;
                            obj.pdf_func{i} = Function('pdf', {x,obj.theta,obj.u}, {p});
                            obj.loglik{i}=log(p);
                            obj.loglik_func{i}=Function('loglik',{x,obj.theta,obj.u});
                            %...FIM and Likelihood incomplete
                        otherwise
                            error('Unexpected distribution type. No PDF created.');
                    end
                end
                
                
            else
                error('Improper number of inputs')
            end
            obj.y_vec = horzcat(obj.y{:});
            
            obj.loglik_tot = 0;
            for i=1:length(obj.loglik)
                obj.loglik_tot=obj.loglik_tot + obj.loglik{i};
            end
            obj.loglik_tot_func = Function('loglik_tot',{obj.y_vec,obj.theta,obj.u},{obj.loglik_tot});
            obj.fim_tot = 0;
            for i=1:length(obj.fim)
                obj.fim_tot=obj.fim_tot + obj.fim{i};
            end
            obj.fim_tot_func = Function('fim_tot',{obj.theta,obj.u},{obj.fim_tot});
            
            
            mu_vec = vertcat(obj.mu{:});
            obj.mu_vec_func = Function('mu_vec',{obj.y_vec,obj.theta,obj.u},{mu_vec});
            pdf_vec = vertcat(obj.pdf{:});
            obj.pdf_vec_func = Function('pdf_vec',{obj.y_vec,obj.theta,obj.u},{pdf_vec});
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Classical Optimal Design Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                    
                    FIM_tot=FIM_tot+xi_k*obj.fim_tot_func(theta_in,u_grid(i));
                    FIM_init=FIM_init+(1/Nu)*full(obj.fim_tot_func(theta_in,u_grid(i)));
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
                    
                    FIM_tot=FIM_tot+obj.fim_tot_func(theta_in,u_k);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Bayesian Optimal Design Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %still needed
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility Functions, for convenience
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function prob = evalMu(varargin)
            import casadi.*
            if nargin == 4
                obj=varargin{1};
                x_in=varargin{2};
                theta_in=varargin{3};
                u_in=varargin{4};
                prob=full(obj.mu_vec_func(theta_in,u_in));
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
                prob=full(obj.pdf_vec_func(x_in,theta_in,u_in));
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
                prob=full(obj.loglik_tot_func(x_in,theta_in,u_in));
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
                prob=full(obj.fim_tot_func(theta_in,u_in));
            else
                error('Improper number of inputs')
            end
        end
        
        
    end
    
end

