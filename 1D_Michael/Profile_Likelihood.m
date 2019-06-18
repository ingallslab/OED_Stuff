addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
Omega=90;
SwarmSize=500;
numExp = 20;
numTrials=200;
u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
inc = length(u)/numExp;
if inc<1
   disp('Too many experiments!');
   return
end 
u = u(1:round(inc):end);    
% [u_opt, w_opt]=Ds_opt([0.5,3,9,3,-18.1,111.7],Omega,u,FIM_comp([0.5,3,9,3,-18.1,111.7],Omega,u,[0.1,0.2]),[0.1,0.2]);
% numSamp=round(w_opt*numTrials*numExp);

% k=ones(1,length(u_opt));
% for j=1:length(u_opt)
%    [~,k(j)]=min(abs(u-u_opt(j)));
% end
% u_opt = u(k);

xVals = cell(1,length(u));
hold on
for i=1:length(u)
    path=strcat('1D_Michael/drive_W90/hist_W=90.000000_u=',num2str(u(i),'%1.6f'),'.txt');
    tmp = load(path,'.txt');
    data = datasample(tmp,200,'Replace',false);
    xVals{1,i} = data;
    scatter(ones(length(xVals{1,i}),1)*u(i),xVals{1,i});
end
hold off
lb = [0.1 0.001 1 1];
ub = [1 5 10 4];
disp('xValues loaded, generating symbolic expressions');
syms = generateOptimSymbols(xVals(1,:),u);
uVals=u;
xvals = vertcat(xVals{:});
lF = Function('lF',{syms.optimVars,syms.fixedparams},{syms.logLik_tot});

flag = 0;
diffe = ub-lb;
x0 = [];
while flag == 0
    r = [diffe(1)*rand(SwarmSize,1)+lb(1) diffe(2)*rand(SwarmSize,1)+lb(2) diffe(3)*rand(SwarmSize,1)+lb(3) diffe(4)*rand(SwarmSize,1)+lb(4)];

    tests =[];
    for i = 1:SwarmSize
        tm = generateOptimVars([r(i,:) -15 100],uVals);
        if tm ~=0
           tests = [tests tm]; 
        end
    end
    if (size(tests,1)~=0)&&(size(tests,2)~=0)&&(~isempty(tests(:,1)))
        a_old = full(lF(tests(:,1),xvals(:)));
        for i = 1:size(tests,2)
            a=full(lF(tests(:,i),xvals(:)));
            index=1;
            if a < a_old && a_old < inf
                index = i;
                flag=1;
            end
            x0 = [x0 tests(:,index)];
        end
    end

end
disp(transpose(x0(end-5:end-2,:)));
disp('multistart analysis complete');
%logLikelihood = Function('logLikelihood', {optimVars,Omega_sym},{logLik_tot});      
err='';
try
    syms.solver.stats();
    solution = syms.solver('x0',x0(:,1),'lbg',syms.lbg,'ubg',syms.ubg,'p',xvals(:));
    mini=full(solution.x(end-5:end-2,1));
    
catch error
    disp(error.identifier);
    err=error.identifier;
end
if strcmp(err,'SWIG:RuntimeError')
    mini=NaN;
end
disp(mini);

figure;
subplot(2,2,1);
hold on
profileLikelihood(1,mini,xVals,uVals);
hold off
% subplot(2,2,2);
% hold on
% profileLikelihood(2,mini,xVals,uVals);
% hold off
% subplot(2,2,3);
% hold on
% profileLikelihood(3,mini,xVals,uVals);
% hold off
% subplot(2,2,4)
% hold on
% profileLikelihood(4,mini,xVals,uVals);
% hold off
disp('Done');

function f = profileLikelihood(Z,mini,xVals,uVals)
    zeta = mini(Z);
    range = (zeta-1):0.01:(zeta+3);
    
    sym = profileSymbols(xVals,uVals,Z);
    out=[];
    kval=[];
    for i=1:length(range)
         m = mini;
         m(Z) = range(i);
         x=vertcat(profileOptimize(xVals, uVals, sym, Z, m));
         if ~isnan(x)
            out = [out x];
            kval=[kval range(i)];
         end
    end
    plot(kval,out);
    f=out;
end

function sym= profileSymbols(xVals,uVals,Z)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi')
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all
    Omega = 60;
    uLow = 0.1;
    uHigh = 0.2;   
    
    x_sym = SX.sym('x_sym');
    u_sym = SX.sym('u_sym');
    Omega_sym = SX.sym('Omega_sym');
    
    a0_sym = SX.sym('a0_sym');
    a_sym = SX.sym('a_sym');
    K_sym = SX.sym('K_sym');
    n_sym = SX.sym('n_sym');
    
    theta_sym=[a0_sym a_sym K_sym n_sym];
    c0_sym = SX.sym('c0_sym');
    c1_sym = SX.sym('c1_sym');
    c_sym=[c0_sym c1_sym];
    par_sym=[theta_sym c_sym];
    
    g=a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)-x_sym;
    g_x=jacobian(g,x_sym);
    
    sigma2=-(a0_sym+a_sym*((u_sym+x_sym).^n_sym)./(K_sym+(u_sym+x_sym).^n_sym)+x_sym)/(2*Omega_sym*g_x);
    pi0=1/(1+exp(-(c0_sym+c1_sym*u_sym)));
    
    x_obs_sym = SX.sym('x_obs_sym');
    phi=(1/(sqrt(2*pi*sigma2)))*exp((-(x_obs_sym-x_sym).^2)./(2*sigma2));
    phi_func = Function('phi_func', {x_obs_sym,x_sym,u_sym,theta_sym,Omega_sym}, {phi});
    
    x_low_sym = SX.sym('x_low_sym');
    x_high_sym = SX.sym('x_high_sym');
    
    Lik=(1-pi0)*phi_func(x_obs_sym,x_low_sym,u_sym,theta_sym,Omega_sym)+pi0*phi_func(x_obs_sym,x_high_sym,u_sym,theta_sym,Omega_sym);
    
    logLik=log(Lik);
    logLik_func = Function('logLik_func', {x_obs_sym,x_low_sym,x_high_sym,u_sym,par_sym,Omega_sym}, {logLik});
    
    logLik_tot = 0;
    
    xvals = cell(size(xVals,2),1);
    for i=1:size(xVals,2)
        for j=1:length(xVals{1,i})
            x_symij = SX.sym(strcat('xvals_',num2str(i),'_',num2str(j)));
            xvals{i}=[xvals{i}; x_symij];
        end
    end
    
    
    xstars_m = []; %mean x values -> monostable region
    constraints=[]; %nonlinear constraint functions
    lbg = []; %nonlinear constraint bounds lower
    ubg = []; %nonlinear constraint bounds upper
    lbw = []; %linear lower bounds
    ubw = []; %linear upper bounds
    
    xstars_h =[]; %mean x values -> upper branch, bistable region
    xstars_l =[]; %mean x values -> lower branch, bistable region
    disp('setting up symbolic expression');
    for i =1:length(uVals)
        u = uVals(i);
        
        if u<uLow
            xstar_i = SX.sym(strcat('xstar_',num2str(i)));
            xstars_m=[xstars_m; xstar_i];

            gstar=a0_sym+a_sym*((u+xstar_i).^n_sym)./(K_sym+(u+xstar_i).^n_sym)-xstar_i;
            constraints = [constraints, {gstar}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:length(xvals{i})
                logLik_tot = logLik_tot + log(phi_func(xvals{i}(j),xstar_i,u,theta_sym,Omega));
            end
        elseif u>uHigh
            xstar_i = SX.sym(strcat('xstar_',num2str(i)));
            xstars_m=[xstars_m; xstar_i];

            gstar=a0_sym+a_sym*((u+xstar_i).^n_sym)./(K_sym+(u+xstar_i).^n_sym)-xstar_i;
            constraints = [constraints, {gstar}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:length(xvals{i})
                logLik_tot = logLik_tot + log(phi_func(xvals{i}(j),xstar_i,u,theta_sym,Omega));
            end
        else
            xstar_h = SX.sym(strcat('xstar_h_',num2str(i)));
            xstars_h=[xstars_h; xstar_h];
            xstar_l = SX.sym(strcat('xstar_l_',num2str(i)));
            xstars_l=[xstars_l; xstar_l];

            constraints = [constraints, {xstar_l - xstar_h}];
            lbg = [lbg; -inf];
            ubg = [ubg; -0.001];

            gstar_h=a0_sym+a_sym*((u+xstar_h).^n_sym)./(K_sym+(u+xstar_h).^n_sym)-xstar_h;
            constraints = [constraints, {gstar_h}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            gstar_l=a0_sym+a_sym*((u+xstar_l).^n_sym)./(K_sym+(u+xstar_l).^n_sym)-xstar_l;
            constraints = [constraints, {gstar_l}];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            for j = 1:length(xvals{i})
                logLik_tot = logLik_tot + logLik_func(xvals{i}(j),xstar_l,xstar_h,u,par_sym,Omega);
            end
        end
    end
    cons = vertcat(constraints{:});
    p = vertcat(xvals{:});
    p = [p; theta_sym(Z)];
    parO = par_sym([1:(Z-1) (Z+1):end]);
    optimVars = [xstars_m; xstars_h; xstars_l; transpose(parO)];
    nlp = struct('x', optimVars, 'f', -logLik_tot, 'g', cons, 'p', p);
    options.error_on_fail = true;
    options.ipopt.max_iter = 200;
%    options.ipopt.fixed_variable_treatment='make_constraint';
%     options.monitor = {'nlp_f','nlp_g'};
    solver = nlpsol('solver','ipopt',nlp,options);
    disp('solver has generated, beginning optimization');
    sym = struct('theta_sym',theta_sym,'optimVars',optimVars, 'logLik_tot',-logLik_tot,'loglikF',Function('logLik_tot',{optimVars,p},{-logLik_tot}),'lbg',lbg,'ubg',ubg,'solver',solver,'fixedparams',p);
end

function min = profileOptimize(xVals, uVals,syms,Z,mini)
    addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
    import casadi.*
    close all     
    xvals = vertcat(xVals{:});
    xvals = [xvals; mini(Z)];
    if ~exist('syms','var')
        syms = profileSymbols(xVals,uVals,Z);
    end
    x0 = profileOptimVars([transpose(mini) -15 100],uVals,Z);
    %disp(transpose(x0(end-5:end-2)));
    %logLikelihood = Function('logLikelihood', {optimVars,Omega_sym},{logLik_tot});      
    err='';
    if x0~=0
        try
            syms.solver.stats();
            solution = syms.solver('x0',x0,'lbg',syms.lbg,'ubg',syms.ubg,'p',xvals);

            min = full(solution.f);
        catch error
            disp(error.identifier);
            err=error.identifier;
        end
        if strcmp(err,'SWIG:RuntimeError')
            min=NaN;
        end
        disp(min);
    else
        min = NaN;
    end
end

function o = profileOptimVars(params,uVals,Z)
    Bounds = [0.1,0.2];
    a0= params(1);
    a = params(2);
    K = params(3);
    n = params(4);
    tol=1e-4;
    g_bnd1_cnt=0;
    maxU=5;
%     while g_bnd1_cnt~=1&&g_bnd1_cnt~=3
%         x_tst=0:tol:maxU;
%         g_bnd1_cnt=sum(abs(diff(sign((a0+a.*((Bounds(1)+x_tst).^n)./(K+(Bounds(1)+x_tst).^n)-x_tst)+eps)./2))>0);
%         maxU=maxU+1;
%     end
%     g_bnd2_cnt=0;
%     maxU=5;
%     while g_bnd2_cnt~=1&&g_bnd2_cnt~=3
%         x_tst=0:tol:maxU;
%         g_bnd2_cnt=sum(abs(diff(sign((a0+a.*((Bounds(2)+x_tst).^n)./(K+(Bounds(2)+x_tst).^n)-x_tst)+eps)./2))>0);
%         maxU=maxU+1;
%     end
%     if (g_bnd1_cnt<3)||(g_bnd2_cnt<3)
%         o=0;
%         return
%     end
    

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
    params(Z)=[];
    o=[x0starM;x0starH;x0starL;transpose(params)];
end