numFits = 50;
numExp = 10;
numTrials = 200;
options = struct('DataSource','Norm_Linspace_LNA','SystemSize',90);
fits = generateFits(numFits,numExp,numTrials,options);

disp(fits);

std_a0 = std(fits(:,1));
std_a = std(fits(:,2));
std_K = std(fits(:,3));
std_n = std(fits(:,4));

disp('Covariance Matrix');
disp(cov(fits));
disp(strcat('Variance in a0: ',num2str((std_a0)^2)));
disp(strcat('Variance in a: ',num2str((std_a)^2)));
disp(strcat('Variance in K: ',num2str((std_K)^2)));
disp(strcat('Variance in n: ',num2str((std_n)^2)));