function min = casadiOptimize(xVals, uVals, lb, ub, SwarmSize, syms)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all     
    xvals = vertcat(xVals(:));
    if ~exist('syms','var')
        syms = generateOptimSymbols(xVals,uVals);
    end
    lF = Function('lF',{syms.optimVars,syms.fixedparams},{syms.logLik_tot});
    
    disp('symbolic expressions generated, now performing multistart analysis');
    
    flag = 0;
    diffe = ub-lb;
    x0 = [];
    while flag == 0
        r = [diffe(1)*rand(SwarmSize,1)+lb(1) diffe(2)*rand(SwarmSize,1)+lb(2) diffe(3)*rand(SwarmSize,1)+lb(3) diffe(4)*rand(SwarmSize,1)+lb(4)];

        tests =[];
        for i = 1:SwarmSize
            tm = generateOptimVars([r(i,:) -15 100],uVals);
            if tm ~=0
               tests = [tests tm]; 
            end
        end
        if size(tests(:,1))~=0
            a_old = full(lF(tests(:,1),xvals));
            for i = 1:size(tests,2)
                a=full(lF(tests(:,i),xvals));
                index=1;
                if a < a_old && a_old < inf
                    a_old = a;
                    index = i;
                    flag = 1;
                end
                x0 = tests(:,index);
            end
        end
        
    end
    disp(transpose(x0(end-5:end-2)));
    disp('multistart analysis complete');
    %logLikelihood = Function('logLikelihood', {optimVars,Omega_sym},{logLik_tot});


    try
        solution = syms.solver('x0',x0,'lbg',syms.lbg,'ubg',syms.ubg,'p',xvals);
        min = full(solution.x(end-5:end-2));
    catch error
        disp(error.identifier);
        if strcmp(error.identifier, 'SWIG:RuntimeError')
           min = NaN; 
        end
    end
    disp(min);
end


function o = generateOptimVars(params,uVals)
    Bounds = [0.1,0.2];
    a0= params(1);
    a = params(2);
    K = params(3);
    n = params(4);
    tol=1e-4;
    g_bnd1_cnt=0;
    maxU=5;
    while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    g_bnd2_cnt=0;
    maxU=5;
    while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
        o=0;
        return
    end
    

    x0starH=[];
    x0starL=[];
    x0starM=[];
    uLow = 0.1;
    uHigh = 0.2;
    theta = params([1,2,3,4]);
    for i=1:length(uVals)
        [lpt,~,hpt] = fixed_point_v4(uVals(i),theta);
        if uVals(i)<uLow
            x0starM=[x0starM;lpt];
        elseif uVals(i)>uHigh
            x0starM=[x0starM;hpt];
        else
            x0starH=[x0starH;hpt];
            x0starL=[x0starL;lpt];
        end
    end

    o=[x0starM;x0starH;x0starL;transpose(params)];
end

