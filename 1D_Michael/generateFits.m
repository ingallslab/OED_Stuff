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
            disp('xValues loaded, generating symbolic expressions');
            syms = generateOptimSymbols(horzcat(xVals{:,1}),u);
            for i=1:numFits
                min = transpose(casadiOptimize(horzcat(xVals{:,i}),u,[0.1 0.001 1 1],[1 5 10 4],100,syms));
                thetas = [thetas; min];
            end
        end
    end
end