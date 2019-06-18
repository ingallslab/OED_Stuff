%   This generates a plot of the loglikelihood function in a neighborhood
%   around the optimal region for a certain dataset. Will help develop an
%   optimized multistart algorithm.
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
[hess,grad] = hessian(syms.logLik_tot,syms.theta_sym);
gradF = Function('gradF',{syms.optimVars,syms.fixedparams},{grad});
hessF = Function('hessF',{syms.optimVars,syms.fixedparams},{hess});
G=full(gradF(full(solution.x(:,1)),xvals));
H=full(hessF(full(solution.x(:,1)),xvals));

subplot(2,2,1);
hold on
plotlik(solution,mini,xvals,syms,1,H);
%plotQuad(solution,mini,xvals,syms,1,G,H);
hold off
subplot(2,2,2);
hold on
plotlik(solution,mini,xvals,syms,2,H);
%plotQuad(solution,mini,xvals,syms,2,G,H);
hold off
subplot(2,2,3);
hold on
plotlik(solution,mini,xvals,syms,3,H);
%plotQuad(solution,mini,xvals,syms,3,G,H);
hold off
subplot(2,2,4);
hold on
plotlik(solution,mini,xvals,syms,4,H);
%plotQuad(solution,mini,xvals,syms,4,G,H);
hold off


function plotlik(solution,mini,xvals,syms,K,H)
    ep = -3:0.01:3;
    [v,~ ]=eig(H);
    
    L =[];
    for i=1:length(ep)
        x0m=full(solution.x(:,1));
        x0m(end-5:end-2)=mini+ep(i).*v(:,K);
        L=[L full(syms.loglikF(x0m,xvals))];
    end
    plot(ep,L);
end

function plotQuad(solution,mini,xvals,syms,K,G,H)

    lower = mini(K)-1; upper= mini(K)+1;
    mVals = lower:0.01:upper;
    L =[];
    x0=full(solution.x(:,1));
    for i=1:length(mVals)
        L=[L quadApprox(x0,mVals(i),mini,K,xvals,syms,G,H)];
    end
    plot(mVals,L);
end

function out=quadApprox(x0,X,mini,i,xvals,syms,G,H)
    x=mini(i);
    x0m=x0;
    x0m(end-6+i)=X;
    out=0;
    out=out+full(syms.loglikF(x0,xvals));
    dx=(x0m-x0);
    out=out+dot(G,dx(end-5:end-2));
    out=out+dot(dx(end-5:end-2),H*dx(end-5:end-2));
end
