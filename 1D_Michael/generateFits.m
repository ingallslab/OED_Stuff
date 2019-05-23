function thetas=generateFits(numFits, numExp, numTrials, options)
    thetas=[];

    if strcmp(options.DataSource,'SSA_Linspace')
        if options.SystemSize == 90
            u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
            u_new=[];
            inc = length(u)/numExp;
            if inc<1
                disp('Too many experiments!');
                return
            end 
            xVals = cell(numExp,numFits);
            for i=1:inc:length(u)
                tmp = load(strcat('drive_W90/hist_W=90.000000_u=',num2str(u(round(i)),'%1.6f'),'.txt'));
 
                for j = 1:numFits
                    data = datasample(tmp,numTrials,'Replace',false);
                    xVals{round(i/inc+1),j}=data;
                end
                if length(u_new)<numExp
                    u_new=[u_new u(round(i))];
                end
            end
            u=u_new;
            disp('xValues loaded');
            for i=1:numFits
                objective = @(theta) -computeLikelihood_v2(cat(2,xVals{:,i}),u,theta);
                options = optimoptions('particleswarm','Display','iter','SwarmSize',200,'FunctionTolerance',1e-3, 'ObjectiveLimit', 0.9e8);
                min = particleswarm(objective,4,[0.001 0.001 1 1],[1 5 10 4],options);
                option2 = optimset('Display','iter');
                min = fminsearch(objective, min,option2);
                thetas = [thetas; min];
            end
        end
        
        
    elseif strcmp(options.DataSource,'SSA_Custom')
    
    else
        
    end  
end