function [low,mid,high]=fixed_point_v4(u,theta)

    a0=theta(1);
    a=theta(2);
    K=theta(3);
    n=theta(4);

    g_func=@(x0,u,theta) theta(1)+theta(2)*((u+x0).^theta(4))./(theta(3)+(u+x0).^theta(4))-x0;
    
    tol=1e-4;
    maxU=5;
    %x_tst=[0:tol:maxU];
    
    sgn_cnt=0;
    while sgn_cnt~=1&&sgn_cnt~=3
        x_tst=[0:tol:maxU];
        gvals=(a0+a.*((u+x_tst).^n)./(K+(u+x_tst).^n)-x_tst)+eps;
        sgn_chng=abs(diff(sign(gvals)./2))>0;
        sgn_cnt=sum(sgn_chng);
        maxU=maxU+1;
    end
    
    x_strt_pts=(x_tst([sgn_chng false]).*gvals([false sgn_chng])...
                -x_tst([false sgn_chng]).*gvals([sgn_chng false]))...
                    ./(gvals([false sgn_chng])-gvals([sgn_chng false]));
    %x_strt_pts=unique(x_strt_pts);
    
    if ~(sgn_cnt==1||sgn_cnt==3)
        low=NaN;
        mid=NaN;
        high=NaN;
        return 
    end

    func=@(x) g_func(x,u,theta);
    low=fzero(func,x_strt_pts(1));
    if(sgn_cnt==1)
        mid=low;
        high=low;
    else
        mid=fzero(func,x_strt_pts(2));
        high=fzero(func,x_strt_pts(3));
    end
end