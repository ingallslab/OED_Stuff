function loglikelihood = computeLikelihood_v2(xVals, uVals, params, s_low,s_high)    
    theta = [params(1), params(2), params(3)];
    c = [params(4), params(5)];
    uLow = 0.1; uHigh=0.2;
    
    
    g = @(u,x,th) theta(1) + (theta(2)*(u+x).^theta(4))/(theta(3)+(u+x).^theta(4)) - x;
    
    [l,~,h]= fixed_point_v4(uLow,theta);
    [L,~,H]= fixed_point_v4(uHigh,theta);
    if (l==h)||(L==H)
        loglikelihood = -1e8;
        return
    end
    
    
    uMean_monostable1=[]; uMean_monostable2=[];uMean_lowbranch=[];uMean_highbranch=[];uMean_midbranch=[];
    xMean_monostable1=[]; xMean_monostable2=[];xMean_lowbranch=[];xMean_highbranch=[];xMean_midbranch=[];
    flg=0;
    for i=1:length(uVals)
        [lPt,mPt,hPt]=fixed_point_v4(uVals(i),theta);

        if (lPt==hPt&&flg==0)
            uMean_monostable1=[uMean_monostable1 experimentInputs(i)];
            xMean_monostable1=[xMean_monostable1 lPt];

        elseif(lPt==hPt&&flg==1)
            uMean_monostable2=[uMean_monostable2 experimentInputs(i)];
            xMean_monostable2=[xMean_monostable2 lPt];
        else
            uMean_lowbranch=[uMean_lowbranch experimentInputs(i)];
            xMean_lowbranch=[xMean_lowbranch lPt];
            uMean_highbranch=[uMean_highbranch experimentInputs(i)];
            xMean_highbranch=[xMean_highbranch hPt];
            uMean_midbranch = [uMean_midbranch experimentInputs(i)];
            xMean_midbranch=[xMean_midbranch mPt];
            flg=1;
        end
    end
    
    likelihood = 0;
    for i=1:length(uVals)
        for j=1:size(xVals,2)
            
        end
    end
    
    
    
end

function ga = gaussLik(mean,sigma,inp)
    ga = (1/(sqrt(2*pi)*sigma))*exp(((mean-inp).^2)./(2*sigma.^2));
end

function out = rho(c,u,B)
    if u <= B(1)
        out = 0;
    elseif u >= B(2)
        out = 1;
    else
        out = 1/(1+exp(-(c(1)+c(2)*u)));
    end
end
