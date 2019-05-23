numFits = 5;
numExp = 10;
numTrials = 20;
options = struct('DataSource','SSA_Linspace','SystemSize',90);
disp(generateFits(numFits,numExp,numTrials,options));