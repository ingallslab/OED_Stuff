function Yend=SSA_Compute_1D(initX,u,OMEGA,finTime,num)
    %initX = initial points
    %u_vals = inputs to simulate
    %OMEGA = system size
    %fintime = how long to simulate for
    %num = number of simulations to do
    if size(initX,1)==1
        X0=repmat(initX,num,1);
    elseif size(initX,1)==num
        X0=initX;
    else
        error('Starting points must be a single initial condition or one for each ensemble instance!')
    end
    
    a0=0.5;
    a=3;
    K=9;
    n=3;

    delta=1.65e-2;

    M = 2; %Number of reactions
    N = 1; %Number of reactants
    
    Yend=zeros(num,1);
    
    R_mu = zeros(M,N);
    %             r1  r2
    R_mu(1) = 1;  % ?b? r1s are made
    R_mu(2) = -1;  % one r1 is lost

    a_mu=zeros(1,M);
    rand('state',sum(100*clock));  %Set the uniform random number generator

    for i=1:num
    
        X=OMEGA*X0(i,:);%zeros(N,5000);t=zeros(1,5000);
        t=0;

        %------------------------ MAIN LOOP ------------------------
        while t<=finTime %main loop

            %step 1: Calculate a_mu & a_0
            % This step calculates the rate at which *any* reaction is expected to occur
            a_mu(1)=delta*OMEGA*(a0 + a*((u+(X(1)/OMEGA))^n)/(K+(u+(X(1)/OMEGA))^n));
            a_mu(2)=delta*X(1);%OMEGA*X(1)/OMEGA;% Degradation of r1
            a_0  = sum(a_mu);


            %Step 2: calculate tau and mu using random number generators
            % Assuming an exponential distribution of reaction rates, the algorithm
            % draws a time for the next reaction to occur (tau) and decides which
            % reaction it will be (next_mu)
            r1  = rand;
            tau = (1/a_0)*log(1/r1);
            next_mu=find(cumsum(a_mu) >= r1*a_0, 1);

            %Step 3: carry out the reaction mu_next
            % Carries out the reaction and advances the time
            t = t + tau;
            Xold=X;
            X=X+ R_mu(next_mu,1:N);
            

        end %end of main loop
        Yend(i)=Xold/OMEGA;
        if mod(i,1000)==0
            disp(i);
        end
    end
    disp('Finished a sim!');
    
end
