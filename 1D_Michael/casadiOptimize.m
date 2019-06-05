function min = casadiOptimize(xVals, uVals, lb, ub, SwarmSize, syms)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all     
    xvals = vertcat(xVals{:});
    if ~exist('syms','var')
        syms = generateOptimSymbols(xVals,uVals);
    end
    lF = Function('lF',{syms.optimVars,syms.fixedparams},{syms.logLik_tot});
       
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
        if (size(tests,1)~=0)&&(size(tests,2)~=0)&&(~isempty(tests(:,1)))
            a_old = full(lF(tests(:,1),xvals(:)));
            for i = 1:size(tests,2)
                a=full(lF(tests(:,i),xvals(:)));
                index=1;
                if a < a_old && a_old < inf
                    index = i;
                    flag=1;
                end
                x0 = [x0 tests(:,index)];
            end
        end
        
    end
    disp(transpose(x0(end-5:end-2,:)));
    disp('multistart analysis complete');
    %logLikelihood = Function('logLikelihood', {optimVars,Omega_sym},{logLik_tot});      
    err='';
    try
        syms.solver.stats();
        solution = syms.solver('x0',x0(:,1),'lbg',syms.lbg,'ubg',syms.ubg,'p',xvals(:));
        
        min = full(solution.x(end-5:end-2,1));
    catch error
        disp(error.identifier);
        err=error.identifier;
    end
    if strcmp(err,'SWIG:RuntimeError')
        min=NaN;
    end

    disp(min);
end




