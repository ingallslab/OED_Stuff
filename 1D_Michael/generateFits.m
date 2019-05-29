function thetas=generateFits(numFits, numExp, numTrials, options)
    thetas=[];
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    Omega = 90;
    
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
    sigma_func = Function('sigma_func',{theta_sym,Omega_sym,u_sym,x_sym},{sqrt(sigma2)});
    
    pi0=1/(1+exp(-(c0_sym+c1_sym*u_sym)));
    pi0_func = Function('pi0_func',{par_sym,u_sym},{pi0});
    
    if strcmp(options.DataSource,'SSA_Linspace')
        if options.SystemSize == 90
            u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
            u_new=[];
            inc = length(u)/numExp;
            if inc<1
                disp('Too many experiments!');
                return
            end 
            xVals = cell(numExp,numFits);
            for i=1:inc:length(u)
                tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u(round(i)),'%1.6f'),'.txt'));
 
                for j = 1:numFits
                    data = datasample(tmp,numTrials,'Replace',false);
                    xVals{round(i/inc+1),j}=data;
                end
                if length(u_new)<numExp
                    u_new=[u_new u(round(i))];
                end
            end
            u=u_new;
            disp('xValues loaded, generating symbolic expressions');
            syms = generateOptimSymbols(horzcat(xVals{:,1}),u);
            lb = [0.1 0.001 1 1];
            ub = [1 5 10 4];
            while numFits > 0
                for i=1:numFits         
                    min = transpose(casadiOptimize(horzcat(xVals{:,i}),u,lb,ub,100,syms));
                    flag=1;
                    for k=1:4
                        if ~(lb(k)<min(k)&&ub(k)>min(k))
                            flag=0;
                            break
                        end
                    end
                    if (~isnan(min))&(flag==1)
                        thetas = [thetas; min];
                    end
                end
                numFits = numFits - size(thetas,1);
            end
        end
    end
    if strcmp(options.DataSource,'Norm_Linspace_LNA')
        if options.SystemSize == 90
            u = linspace(0,0.3,numExp);
            xVals=cell(numExp,numFits);
            theta_t=[0.5,3,9,3];
            for i=1:numExp
                [lpt,~,hpt]=fixed_point_v4(u(i),theta_t);
                
                flag = 0;
                if flag==0&&lpt==hpt
                    for j=1:numFits 
                        sigma=full(sigma_func(theta_t,Omega,u(i),lpt));
                        xVals{i,j}=normrnd(lpt,sigma,[numTrials,1]);   
                    end
                elseif flag==1&&lpt==hpt
                    for j=1:numFits 
                        sigma=full(sigma_func(theta_t,Omega,u(i),hpt));
                        xVals{i,j}=normrnd(hpt,sigma,[numTrials,1]);   
                    end
                else
                    flag=1;
                    for j=1:numFits 
                        sigmah=full(sigma_func(theta_t,Omega,u(i),hpt));
                        sigmal=full(sigma_func(theta_t,Omega,u(i),lpt));
                        numHigh = round(full(pi0_func([-18.1,111.7],u(i)))*numTrials);
                        numLow = numTrials-numHigh;
                        xVals{i,j}=[normrnd(hpt,sigmah,[numHigh,1]); normrnd(lpt,sigmal,[numLow,1])];
                    end
                end
            end

            
            disp('xValues loaded, generating symbolic expressions');
            syms = generateOptimSymbols(horzcat(xVals{:,1}),u);
            lb = [0.1 0.001 1 1];
            ub = [1 5 10 4];
            while numFits > 0
                for i=1:numFits
                    min = transpose(casadiOptimize(horzcat(xVals{:,i}),u,lb,ub,500,syms));
                    flag = 1;
                    for k=1:4
                        if ~(lb(k)<min(k)&&ub(k)>min(k))
                            flag=0;
                            break
                        end
                    end
                    if (~isnan(min))&(flag==1)
                        thetas = [thetas; min];
                    end
                end
                numFits = numFits - size(thetas,1);
            end
        end
    end
end