function o = generateOptimVars(params,uVals)
    Bounds = [0.1,0.2];
    a0= params(1);
    a = params(2);
    K = params(3);
    n = params(4);
    tol=1e-4;
    g_bnd1_cnt=0;
    maxU=5;
    while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    g_bnd2_cnt=0;
    maxU=5;
    while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
        x_tst=0:tol:maxU;
        g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
        maxU=maxU+1;
    end
    if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
        o=0;
        return
    end
    

    x0starH=[];
    x0starL=[];
    x0starM=[];
    uLow = 0.1;
    uHigh = 0.2;
    theta = params([1,2,3,4]);
    for i=1:length(uVals)
        [lpt,~,hpt] = fixed_point_v4(uVals(i),theta);
        if uVals(i)<uLow
            x0starM=[x0starM;lpt];
        elseif uVals(i)>uHigh
            x0starM=[x0starM;hpt];
        else
            x0starH=[x0starH;hpt];
            x0starL=[x0starL;lpt];
        end
    end

    o=[x0starM;x0starH;x0starL;transpose(params)];
end