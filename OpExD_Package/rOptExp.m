classdef rOptExp
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xi_opt
        u_opt    
    end
    
    methods
        function obj=rOptExp(u,xi)
            obj.xi_opt=xi;
            obj.u_opt=u;
        end
        
        function eoe=genExact(varargin)
            obj=varargin{1};
            N=varargin{2};
            u_ext=[];
            xiSum=0;
            repSum=0;
            for i=1:length(obj.xi_opt)
                xiSum=xiSum+N*obj.xi_opt(i);
                reps=round(xiSum-repSum);
                if round(xiSum-repSum)>0
                    u_ext=[u_ext repmat(obj.u_opt(i),1,reps)];
                end
                repSum=repSum+reps;
            end
            
            eoe=eOptExp(u_ext);
        end
        
    end
    
end

