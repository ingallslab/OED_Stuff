function loglikelihood = computeLikelihood(xVals, uVals, theta, sigmaVals)
    loglikelihood = 0;
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    %casadi setup
    x_sym = SX.sym('x');
    u_sym = SX.sym('u');
    a_sym = SX.sym('a');
    K_sym = SX.sym('K');
    n_sym = SX.sym('n');
    theta_sym = [a_sym,K_sym,n_sym];
    sigma_sym = SX.sym('sigma');

    xstar = a_sym *(u_sym.^n_sym)./(u_sym.^n_sym+K_sym.^n_sym); 

    x_obs_sym = SX.sym('x_obs');

    lik = (1/(sqrt(2*pi)*sigma_sym))*exp(-((x_obs_sym-xstar).^2)./(2*sigma_sym.^2));
    loglik = log(lik);
    loglik_f = Function('loglik_f',{x_obs_sym,u_sym,theta_sym,sigma_sym},{loglik});
    
    %dataMatrix = [uVals;sigmaVals;xVals];
    %disp(dataMatrix);
    for i=1:size(xVals,1)
        for j=1:size(xVals,2)
            if sigmaVals(j)~=0
                temp = loglik_f(xVals(i,j),uVals(j),theta,sigmaVals(j));
                loglikelihood = loglikelihood + temp;
            end
            
            %disp(loglikelihood);
        end
    end
    
    loglikelihood=full(loglikelihood);
   
end