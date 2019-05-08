addpath('/Users/mrastwoo/Documents/MATLAB/Casadi/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
close all

%some nice values
a=3;
K=9;
n=3;
m=2;

%Casadi Setup
%%variables
x_sym = SX.sym('x');
u_sym = SX.sym('u');
a_sym = SX.sym('a');
K_sym = SX.sym('K');
n_sym = SX.sym('n');
theta = [a_sym,K_sym,n_sym];

g = a_sym*(u_sym.^n_sym)./(K_sym.^n_sym+u_sym.^n_sym)-x_sym;

%%plotting steady states -> plots x* as function of u starting with casadi
%%%%%%%%%%%%%%%%%%%%%%%%%%%

xstar = a_sym *(u_sym.^n_sym)./(u_sym.^n_sym+K_sym.^n_sym); 
%casadi variable for steady state

xStar = Function('xStar',{u_sym,theta},{xstar}); 
%tell matlab to consider xStar as a function with value xstar

uVals = linspace(0,40,100);
sVals = xStar(uVals,[a,K,n]); %s stands for steady state
zVals = xStar(uVals,[a,K,m]); %I wanted to see what n=2 looks like

%sparse converts casadi.DM classtype to a list of doubles
plot(uVals,sparse(sVals)); 
hold on
plot(uVals,sparse(zVals));
hold off

%generating normally distributed data points (~N(0,0.1))
disp(length(sparse(sVals)));
normRands = normrnd(0,0.1,[100,1]);
Redist_sVals = abs(normRands+sVals(:));%normally redistributed sVals
scatter(uVals,sparse(Redist_sVals),'filled');

%generating more interesting normally distributed data points
xNormRands = linspace(0,1,100); %preallocating array for efficiency
%x dependent random variables -> N(0,x*0.2)
for i=1:100
    xNormRands(i) = normrnd(0,0.2*sVals(i));
end
xRedist_sVals = abs(xNormRands+sVals(:));
scatter(uVals,sparse(xRedist_sVals),'filled');


%generating data points from SSA -> more complicated!