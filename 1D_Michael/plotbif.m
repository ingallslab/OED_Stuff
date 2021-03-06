function plotbif(uVals, theta, options)
    u_monostable1=[]; u_monostable2=[];u_lowbranch=[];u_highbranch=[];u_midbranch=[];
    x_monostable1=[]; x_monostable2=[];x_lowbranch=[];x_highbranch=[];x_midbranch=[];
    flg=0;
    for i=1:length(uVals)
        [lPt,mPt,hPt]=fixed_point_v4(uVals(i),theta);

        if (lPt==hPt&&flg==0)
            u_monostable1=[u_monostable1 uVals(i)];
            x_monostable1=[x_monostable1 lPt];

        elseif(lPt==hPt&&flg==1)
            u_monostable2=[u_monostable2 uVals(i)];
            x_monostable2=[x_monostable2 lPt];
        else
            u_lowbranch=[u_lowbranch uVals(i)];
            x_lowbranch=[x_lowbranch lPt];
            u_highbranch=[u_highbranch uVals(i)];
            x_highbranch=[x_highbranch hPt];
            u_midbranch = [u_midbranch uVals(i)];
            x_midbranch=[x_midbranch mPt];
            flg=1;
        end
    end
    U_Low=[u_monostable1 u_lowbranch];
    X_Low=[x_monostable1 x_lowbranch];
    U_High=[u_highbranch u_monostable2];
    X_High=[x_highbranch x_monostable2];
    
    
        
    if exist('options','var')&&strcmp(options, 'black')
        plot(U_Low,X_Low,'k');
        plot(u_midbranch,x_midbranch,'k:');
        plot(U_High,X_High,'k');
    else
        plot(U_Low,X_Low,'b','LineWidth',2);
        plot(u_midbranch,x_midbranch,'k:','LineWidth',2);
        plot(U_High,X_High,'r','LineWidth',2);
    end
    
end

