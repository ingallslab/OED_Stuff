function loglikelihood = computeLikelihood_v2(xVals, uVals, theta,sigma)    
    uLow = 0.11; uHigh=0.19;
    
    a0= theta(1);
    a = theta(2);
    K = theta(3);
    n = theta(4);
    Omega = 90;
    tol=1e-4;
    maxU=5;
    for i =1:length(uVals)
        sgn_cnt=0;
        u = uVals(i);
        while sgn_cnt~=1&&sgn_cnt~=3
            x_tst=[0:tol:maxU];
            gvals=(a0+a.*((u+x_tst).^n)./(K+(u+x_tst).^n)-x_tst)+eps;
            sgn_chng=abs(diff(sign(gvals)./2))>0;
            sgn_cnt=sum(sgn_chng);
            maxU=maxU+1;
        end
        if (sgn_cnt ~= 3)&&(uVals(i) <= uHigh)&&(uVals(i) >= uLow)
            loglikelihood = -1e8;
            return
        end
    end
    
    if exist('sigma','var')
        format long;
        likelihood = 0;
        for i=1:length(uVals)
            u = uVals(i);
            [lpt,~,hpt]=fixed_point_v4(u,theta);
            temp=0;
            if(sigma(1,i)~=0)&&(sigma(2,i)==0)
                temp=sum(log(...
                normpdf(xVals(:,i),lpt,sigma(1,i)) ));
                likelihood = likelihood +temp;
            elseif (sigma(1,i)~=0)&&(sigma(2,i)~=0)
    %             disp([sigma(1,i) sigma(2,i)]);
    %             disp([lpt hpt]);
                temp=(1-rho(u))*normpdf(xVals(:,i),lpt,sigma(1,i))+...
                rho(u)*normpdf(xVals(:,i),hpt,sigma(2,i)) ;
    %             disp([normpd(xVals(:,i),lpt,sigma(1,i)) normpd(xVals(:,i),hpt,sigma(2,i)) xVals(:,i)]);
                for j = 1:length(xVals(:,i))
                    if temp(j)~=0
                        likelihood=likelihood+log(temp(j));
                    end
                end
            elseif(sigma(2,i)~=0)&&(sigma(1,i)==0)
                temp=sum(log(...
                normpdf(xVals(:,i),hpt,sigma(2,i)) ));
                likelihood = likelihood +temp;
            else
                likelihood=-1e9;
            end

        end   
        if likelihood==inf
           likelihood = -1e8; 
        end
        loglikelihood = likelihood;
    else
        format long;
        likelihood = 0;
        for i=1:length(uVals)
            u = uVals(i);
            [lpt,~,hpt]=fixed_point_v4(u,theta);
            if u < uLow
                sigmalow=(a0+a*((u+lpt)^n)/(K+(u+lpt)^n)+lpt)/sqrt(Omega);
                temp=sum(log(...
                normpdf(xVals(:,i),lpt,sigmalow) ));
                likelihood = likelihood +temp;
            elseif uLow <= u && u <= uHigh
                sigmalow=(a0+a*((u+lpt)^n)/(K+(u+lpt)^n)+lpt)/sqrt(Omega);
                sigmahigh=(a0+a*((u+hpt)^n)/(K+(u+hpt)^n)+hpt)/sqrt(Omega);
                temp=(1-rho(u))*normpdf(xVals(:,i),lpt,sigmalow)+...
                rho(u)*normpdf(xVals(:,i),hpt,sigmahigh) ;
                for j = 1:length(xVals(:,i))
                    if temp(j)~=0
                        likelihood=likelihood+log(temp(j));
                    end
                end
            elseif uHigh < u
                sigmahigh=(a0+a*((u+hpt)^n)/(K+(u+hpt)^n)+hpt)/sqrt(Omega);
                temp=sum(log(...
                normpdf(xVals(:,i),hpt,sigmahigh) ));
                likelihood = likelihood +temp;
            else
                likelihood=-1e9;
            end

        end   
        if likelihood==inf
           likelihood = -1e8; 
        end
        loglikelihood = likelihood;
        
    end
end

function out = rho(u)
    B = [0.1 0.2];
    c = [-15,100];
    if u < B(1)
        out = 0;
    elseif u > B(2)
        out = 1;
    else
        out = 1/(1+exp(-(c(1)+c(2)*u)));
    end
end
