function [low,mid,high]=fixed_point_v5(func,maxX)
  
    tol=1e-2;
    
    increment=round(.2*maxX);
    
    sgn_cnt=0;
    while sgn_cnt~=1&&sgn_cnt~=3
        x_tst=[0:tol:maxX];
        gvals=func(x_tst)+eps;
        sgn_chng=abs(diff(sign(gvals)./2))>0;
        sgn_cnt=sum(sgn_chng);
        maxX=maxX+increment;
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

    low=fzero(func,x_strt_pts(1));
    if(sgn_cnt==1)
        mid=low;
        high=low;
    else
        mid=fzero(func,x_strt_pts(2));
        high=fzero(func,x_strt_pts(3));
    end
end