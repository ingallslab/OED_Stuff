function [nans,thetas]=generateFits(numFits, numExp, numTrials, options)
    thetas=[];
    
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    Omega = options.SystemSize;
    par_true = [0.5,3,9,3,-18.1,111.7];
    
    
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
    numNans=-0;
    if strcmp(options.DataSource,'SSA')
        if options.SystemSize==60
            loadData = dlmread('SSA_Data_60.txt');
            if strcmp(options.Optimality,'D_Optimal')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:round(inc):end);    
                [u_opt, w_opt]=D_opt_c(par_true,Omega,u,FIM_comp(par_true,Omega,u,[0.1,0.2]),[0.1,0.2]);
                numSamp=round(w_opt*numTrials*numExp);
                
                k=ones(1,length(u_opt));
                for j=1:length(u_opt)
                    [~,k(j)]=min(abs(u-u_opt(j)));
                end
                u_opt=u(k);
                xVals = cell(numFits,length(u_opt));
                
                for i=1:length(u_opt)
                    for j = 1:numFits
                        data = datasample(loadData(:,i),numSamp(i),'Replace',false);
                        xVals{j,i} = data;
                    end
                end
                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals(1,:),u_opt);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits >= 0
                    for i=1:numFits         
                        mini = transpose(casadiOptimize(xVals(i,:),u_opt,lb,ub,100,syms));
                        flag=1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
            elseif strcmp(options.Optimality,'Ds_Optimal')
                
            elseif strcmp(options.Optimality,'LinSpace')
                
            end
        elseif options.SystemSize == 90
            if strcmp(options.Optimality,'D_Optimal')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:round(inc):end);    
                [u_opt, w_opt]=D_opt_c(par_true,Omega,u,FIM_comp(par_true,Omega,u,[0.1,0.2]),[0.1,0.2]);
                numSamp=round(w_opt*numTrials*numExp);
                
                k=ones(1,length(u_opt));
                for j=1:length(u_opt)
                    [~,k(j)]=min(abs(u-u_opt(j)));
                end
                u_opt=u(k);
                xVals = cell(numFits,length(u_opt));
                for i=1:length(u_opt)
                    tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u_opt(i),'%1.6f'),'.txt'));
                    for j = 1:numFits
                        data = datasample(tmp,numSamp(i),'Replace',false);
                        xVals{j,i} = data;
                    end
                end
                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals(1,:),u_opt);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits >= 0
                    for i=1:numFits         
                        mini = transpose(casadiOptimize(xVals(i,:),u_opt,lb,ub,100,syms));
                        flag=1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
                
                
            elseif strcmp(options.Optimality,'Ds_Optimal')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:inc:end);        
                [u_opt,w_opt] = Ds_opt(par_true,Omega,u,FIM_comp(par_true,Omega,u,[0.1,0.2]),[0.1,0.2]);
                numSamp=round(w_opt*numTrials*numExp);
                
                k=ones(1,length(u_opt));
                for j=1:length(u_opt)
                    [~,k(j)]=min(abs(u-u_opt(j)));
                end
                u_opt=u(k);
                
                xVals = cell(numFits,length(u_opt));
                for i=1:length(u_opt)
                    tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u_opt(i),'%1.6f'),'.txt'));
                    for j = 1:numFits
                        data = datasample(tmp,numSamp(i),'Replace',false);
                        xVals{j,i} = data;
                    end
                end
                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals(1,:),u_opt);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits >= 0
                    for i=1:numFits         
                        mini = transpose(casadiOptimize(xVals(i,:),u_opt,lb,ub,100,syms));
                        flag=1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
                
                    
                
            elseif strcmp(options.Optimality, 'LinSpace')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                u_new=[];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                xVals = cell(numFits,numExp);
                for i=1:inc:length(u)
                    tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u(round(i)),'%1.6f'),'.txt'));

                    for j = 1:numFits
                        data = datasample(tmp,numTrials,'Replace',false);
                        xVals{j,round(i/inc)+1} = data;
                    end
                    if length(u_new)<numExp
                        u_new=[u_new u(round(i))];
                    end
                end
                u=u_new;
                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals(1,:),u);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits >= 0
                    for i=1:numFits         
                        mini = transpose(casadiOptimize(xVals(i,:),u,lb,ub,100,syms));
                        flag=1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
            end
        end
    end
    if strcmp(options.DataSource,'Normal')
        if options.SystemSize == 90
            if strcmp(options.Optimality,'D_Optimal')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:inc:end);        
                [u_opt,w_opt] = D_opt_c(par_true,Omega,u,FIM_comp(par_true,Omega,u,[0.1,0.2]),[0.1,0.2]);
                numSamp=round(w_opt*numTrials*numExp);
                
                k=ones(1,length(u_opt));
                for j=1:length(u_opt)
                    [~,k(j)]=min(abs(u-u_opt(j)));
                end
                u_opt=u(k);
                
                xVals = cell(numFits,length(u_opt));
                theta_t=[0.5,3,9,3];
                for i=1:length(u_opt)
                    [lpt,~,hpt]=fixed_point_v4(u_opt(i),theta_t);

                    flag = 0;
                    if flag==0&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u_opt(i),lpt));
                            xVals{j,round(i)}=normrnd(lpt,sigma,[numSamp(i),1]);   
                        end
                    elseif flag==1&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u_opt(i),hpt));
                            xVals{j,round(i)}=normrnd(hpt,sigma,[numSamp(i),1]);   
                        end
                    else
                        flag=1;
                        for j=1:numFits 
                            sigmah=full(sigma_func(theta_t,Omega,u_opt(i),hpt));
                            sigmal=full(sigma_func(theta_t,Omega,u_opt(i),lpt));
                            numHigh = round(full(pi0_func([-18.1,111.7],u_opt(i)))*numSamp(i));
                            numLow = numSamp(i)-numHigh;
                            xVals{j,round(i)}=[normrnd(hpt,sigmah,[numHigh,1]); normrnd(lpt,sigmal,[numLow,1])];
                        end
                    end
                end


                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals,u_opt);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits > 0
                    for i=1:numFits
                        mini = transpose(casadiOptimize(xVals(i,:),u_opt,lb,ub,500,syms));
                        flag = 1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
            elseif strcmp(options.Optimality,'Ds_Optimal')
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:inc:end);        
                [u_opt,w_opt] = Ds_opt(par_true,Omega,u,FIM_comp(par_true,Omega,u,[0.1,0.2]),[0.1,0.2]);
                numSamp=round(w_opt*numTrials*numExp);
                
                k=ones(1,length(u_opt));
                for j=1:length(u_opt)
                    [~,k(j)]=min(abs(u-u_opt(j)));
                end
                u_opt=u(k);
                
                xVals = cell(numFits,length(u_opt));
                theta_t=[0.5,3,9,3];
                for i=1:length(u_opt)
                    [lpt,~,hpt]=fixed_point_v4(u_opt(i),theta_t);

                    flag = 0;
                    if flag==0&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u_opt(i),lpt));
                            xVals{j,round(i)}=normrnd(lpt,sigma,[numSamp(i),1]);   
                        end
                    elseif flag==1&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u_opt(i),hpt));
                            xVals{j,round(i)}=normrnd(hpt,sigma,[numSamp(i),1]);   
                        end
                    else
                        flag=1;
                        for j=1:numFits 
                            sigmah=full(sigma_func(theta_t,Omega,u_opt(i),hpt));
                            sigmal=full(sigma_func(theta_t,Omega,u_opt(i),lpt));
                            numHigh = round(full(pi0_func([-18.1,111.7],u_opt(i)))*numSamp(i));
                            numLow = numSamp(i)-numHigh;
                            xVals{j,round(i)}=[normrnd(hpt,sigmah,[numHigh,1]); normrnd(lpt,sigmal,[numLow,1])];
                        end
                    end
                end


                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals,u_opt);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits > 0
                    for i=1:numFits
                        mini = transpose(casadiOptimize(xVals(i,:),u_opt,lb,ub,500,syms));
                        flag = 1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end    
            elseif strcmp(options.Optimality,'LinSpace')               
                u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
                inc = length(u)/numExp;
                if inc<1
                    disp('Too many experiments!');
                    return
                end 
                u = u(1:inc:end);
                xVals=cell(numFits,numExp);
                theta_t=[0.5,3,9,3];
                for i=1:numExp
                    [lpt,~,hpt]=fixed_point_v4(u(round(i)),theta_t);

                    flag = 0;
                    if flag==0&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u(round(i)),lpt));
                            xVals{j,round(i)}=normrnd(lpt,sigma,[numTrials,1]);   
                        end
                    elseif flag==1&&lpt==hpt
                        for j=1:numFits 
                            sigma=full(sigma_func(theta_t,Omega,u(round(i)),hpt));
                            xVals{j,round(i)}=normrnd(hpt,sigma,[numTrials,1]);   
                        end
                    else
                        flag=1;
                        for j=1:numFits 
                            sigmah=full(sigma_func(theta_t,Omega,u(round(i)),hpt));
                            sigmal=full(sigma_func(theta_t,Omega,u(round(i)),lpt));
                            numHigh = round(full(pi0_func([-18.1,111.7],u(round(i))))*numTrials);
                            numLow = numTrials-numHigh;
                            xVals{j,round(i)}=[normrnd(hpt,sigmah,[numHigh,1]); normrnd(lpt,sigmal,[numLow,1])];
                        end
                    end
                end


                disp('xValues loaded, generating symbolic expressions');
                syms = generateOptimSymbols(xVals,u);
                lb = [0.1 0.001 1 1];
                ub = [1 3.75 20 3.75];
                while numFits > 0
                    for i=1:numFits
                        mini = transpose(casadiOptimize(xVals(i,:),u,lb,ub,500,syms));
                        flag = 1;
                        for k=1:4
                            if ~(lb(k)<mini(k)&&ub(k)>mini(k))
                                flag=0;
                                break
                            end
                        end
                        if (~isnan(mini))&(flag==1)
                            thetas = [thetas; mini];
                            disp(i);
                        else
                            numNans=numNans+1;
                        end
                    end
                    numFits = numFits - size(thetas,1);
                end
            end
        end
    end
    disp(strcat('numNans:',num2str(numNans)));
    nans=numNans;
end