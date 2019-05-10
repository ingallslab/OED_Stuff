function loglikelihood = compLikelihood_nocasadi(xVals, uVals, theta, sigmaVals)
    loglikelihood = 0;
    %dataMatrix = [uVals;sigmaVals;xVals];
    %disp(dataMatrix);
    for i=1:size(xVals,1)
        for j=1:size(xVals,2)
            if sigmaVals(j)~=0
                loglikelihood = loglikelihood + indLik(xVals(i,j),uVals(j),theta,sigmaVals(j));
                
            end
        end
    end
end

function xstar = xStar(theta,u)
    a = theta(1);
    K = theta(2);
    n = theta(3);
    xstar = a*(u.^n)./(u.^n + K.^n);
end


function indlik = indLik(xVal,uVal,theta,sigma)
    indlik=-log(sqrt(2*pi)*sigma)-(((xVal-xStar(theta,uVal)).^2)./(2*sigma.^2));
end