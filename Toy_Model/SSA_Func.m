function YEnd = SSA_Func(initX,u,OMEGA,finTime,num)
    %simulates data from stochastic model by solving master equation
    %iteratively
    if size(initX,1)==1
        X0=repmat(initX,num,1);
    elseif size(initX,1)==num
        X0=initX;
    else
        error('Starting points must be a single initial condition or one for each ensemble instance!')
    end
    a=3;
    K=9;
    n=3;
    delta = 1.65e-2; %approximately 1/60
    
    YEnd = zeros(num,1); %initialize output matrix
    
    R_mu = [1, -1]; %reaction end result possibilities(stoichiometry)
    a_mu = [0 0]; %reaction propensities (times dt)
    
    rand('state',sum(100*clock)); %RNG with current datetime as seed
    for i=1:num
        X = OMEGA*X0(i);%number of molecules generated
        t=0;%time
        
        while t<=finTime
            %recall reaction propensity (transition probabilities): 
            %   Production -> Omega(au^n/(K^n+u^n))
            %   Decay      -> X
            %this generates the probability that the system evolves into
            %each possible next state
            a_mu(1) = (delta*OMEGA*a*u.^n)./(u.^n+K.^n);
            a_mu(2) = delta*X;
            a_0 = sum(a_mu);
            
            %recall that the exponential distribution generates
            %probabilities of time between events in a poisson process
            r = rand; 
            %   generate a random number
            tau = log(1/r)/a_0;
            %   exponential distribution
            next_mu = find(cumsum(a_mu)>=r*a_0,1); 
            %   nice way to pick a random reaction from the random number
            
            t=t+tau; %time tends to increase
            tempX=X;
            X=X+R_mu(next_mu);
        end
        YEnd(i) = tempX/OMEGA;
    end
end

