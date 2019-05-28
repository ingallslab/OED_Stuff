% numFits = 10;
% numExp = 20;
% numTrials = 200;
% options = struct('DataSource','SSA_Linspace','SystemSize',90); %SSA_Linspace, Norm_Linspace_LNA
% fits = generateFits(numFits,numExp,numTrials,options);
% 
% disp(fits);
% 
% std_a0 = std(fits(:,1));
% std_a = std(fits(:,2));
% std_K = std(fits(:,3));
% std_n = std(fits(:,4));
% a_a0 = mean(fits(:,1));
% a_a = mean(fits(:,2));
% a_K = mean(fits(:,3));
% a_n = mean(fits(:,4));
% 
% 
% disp('Covariance Matrix');
% disp(cov(fits));
% disp(strcat("Variance in a0: ",num2str((std_a0)^2)));
% disp(strcat("Variance in a: ",num2str((std_a)^2)));
% disp(strcat("Variance in K: ",num2str((std_K)^2)));
% disp(strcat("Variance in n: ",num2str((std_n)^2)));
% 
% disp(strcat("Mean a0: ",num2str(a_a0)));
% disp(strcat("Mean a: ",num2str(a_a)));
% disp(strcat("Mean K: ",num2str(a_K)));
% disp(strcat("Mean n: ",num2str(a_n)));
