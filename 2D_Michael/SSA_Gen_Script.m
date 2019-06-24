umax_1=1;
umax_2=100;
[u1_grid,u2_grid] = meshgrid(0:0.1:umax_1,0:10:umax_2);
OMEGA=90;
finTime = 5000;
num = 20;
SSAData=cell(size(u1_grid,1),size(u1_grid,2));
for i=1:size(u1_grid,1)
    disp(strcat('i = ',num2str(i)));
    parfor j=1:size(u2_grid,2)
        SSAData{i,j} = Toggle_Ensemble([0.1,0.1],[u1_grid(i,j) u2_grid(i,j)],OMEGA,finTime,num);
        dlmwrite(strcat('2D_Michael/Data/SSAData_',num2str(i),'_',num2str(j)),SSAData{i,j},'\t');
    end
end
function Yend=Toggle_Ensemble(initX,u_vals,OMEGA,finTime,num)
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

    theta=[alpha_1 alpha_2 beta_1 beta_2 K_1 K_2 n_1 n_2 kappa_1 kappa_2 m_1 m_2];

    M = 4; %Number of reactions
    N = 2; %Number of reactants
    
    Yend=zeros(num,N);
    
    R_mu = zeros(M,N);
    %             r1  r2
    R_mu(1,:) = [ 1   0] ;  % ?b? r1s are made
    R_mu(2,:) = [ 0   1] ;  % ?b? r2s are made
    R_mu(3,:) = [ -1  0] ;  % one r1 is lost
    R_mu(4,:) = [ 0  -1] ;  % one r2 is lost

    a_mu=zeros(1,M);
    rand('state',sum(100*clock));  %Set the uniform random number generator

    for i=1:num
        disp(i);
        X=OMEGA*X0(i,:);%zeros(N,5000);t=zeros(1,5000);
        t=0;

        %------------------------ MAIN LOOP ------------------------
        k = 2;
        while t<=finTime %main loop

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
            t = t + tau;
            Xold=X;
            X=X+ R_mu(next_mu,1:N);

        end %end of main loop
        
        Yend(i,:)=Xold/OMEGA;
    end
    
end
