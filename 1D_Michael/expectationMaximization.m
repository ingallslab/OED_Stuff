function estimate = expectationMaximization(xVals,uVals,theta_i)
    uLow = 0.1; uHigh = 0.2;
    for i=1:length(uVals)
        u = uVals(i);
        [lpt,~,hpt]=fixed_point_v4(u,theta_i);
        temp=0;
        numComponents = 1;
        if u >= uLow && u <= uHigh
            numComponents = 2;
        end
        R=[];%responsibility matrix
        for n=1:numComponents
            for j=1:size(xVals,2)
                
            end
        end
            

    end   
    loglikelihood = likelihood;
end