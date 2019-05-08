function t=Toggle_SSA_SwitchTime(initX,u_vals,OMEGA,sddl_pt,ortho_vec,maxTime)

    alpha_1 = 13.609;
    alpha_2 = 60.882;   
    beta_1 = 3529.923;   
    beta_2 = 1053.916;      
    K_1 = 31.94;    
    K_2 = 30.0;    
    n_1  = 2.00;     
    n_2 = 2.00;     
    kappa_1 = 0.0906;   
    kappa_2 = 11.65;    
    m_1= 2.00;     
    m_2 = 2.00;

    delta=1.65e-2;

    X=OMEGA*initX;%zeros(N,5000);t=zeros(1,5000);
    
    dir=sign(initX*ortho_vec'-sddl_pt*ortho_vec');

    M = 4; %Number of reactions
    N = 2; %Number of reactants

    R_mu = zeros(M,N);
    R_mu(1,:) = [ 1   0] ;  % ?b? r1s are made
    R_mu(2,:) = [ 0   1] ;  % ?b? r2s are made
    R_mu(3,:) = [ -1  0] ;  % one r1 is lost
    R_mu(4,:) = [ 0  -1] ;  % one r2 is lost

    a_mu=zeros(1,M);
    %[R_mu] = defineReactions(N,M) ;
    rand('state',sum(100*clock));  %Set the uniform random number generator

    t=0;

    k = 2;
    while and((sign((X/OMEGA)*ortho_vec'-sddl_pt*ortho_vec')==dir), t<maxTime) %main loop
        %(X/OMEGA)*ortho_vec'-sddl_pt*ortho_vec'
        %step 1: Calculate a_mu & a_0
        % This step calculates the rate at which *any* reaction is expected to occur

        %a_mu(1)=OMEGA*am*(1+(X(2)/OMEGA/Kr/2)/w*(2+(X(2)/OMEGA/Kr/2)))/(1+(X(2)/OMEGA/Kr/2))^2; % Synthesis of mRNA for r1
        a_mu(1)=delta*OMEGA*(alpha_1 + beta_1./(1+(((X(2)/OMEGA)./K_2)*(1./(1+(u_vals(2)./kappa_2).^m_2))).^n_1));

        %a_mu(2)=OMEGA*am*(1+(X(1)/OMEGA/Kr/2)/w*(2+(X(1)/OMEGA/Kr/2)))/(1+(X(1)/OMEGA/Kr/2))^2; % Synthesis of mRNA for r2
        a_mu(2)=delta*OMEGA*(alpha_2 + beta_2./(1+(((X(1)/OMEGA)./K_1)*(1./(1+(u_vals(1)./kappa_1).^m_1))).^n_2));
        a_mu(3)=delta*X(1);%OMEGA*X(1)/OMEGA;% Degradation of r1
        a_mu(4)=delta*X(2);%OMEGA*X(2)/OMEGA;% Degradation of r2
        a_0  = sum(a_mu);

        %Step 2: calculate tau and mu using random number generators
        % Assuming an exponential distribution of reaction rates, the algorithm
        % draws a time for the next reaction to occur (tau) and decides which
        % reaction it will be (next_mu)
        r1  = rand;
        tau = (1/a_0)*log(1/r1);
        r2 = rand;
        next_mu=find(cumsum(a_mu) >= r2*a_0, 1);

        %Step 3: carry out the reaction mu_next
        % Carries out the reaction and advances the time
        t =t + tau;
        X=X+ R_mu(next_mu,1:N);

    end; %end of main loop
    
    if t>maxTime
        t=NaN;
    end
    
end