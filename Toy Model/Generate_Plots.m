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

xstar = a_sym *(u_sym.^n_sym)./(u_sym.^n_sym+K_sym.^n_sym); 
%casadi variable for steady state

xStar = Function('xStar',{u_sym,theta},{xstar}); 
%tell matlab to consider xStar as a function with value xstar



%%
%%%-----------------plot generation begins here-----------------
%%%
uVals = linspace(0,40,100);
sVals = xStar(uVals,[a,K,n]); %s stands for steady state
zVals = xStar(uVals,[a,K,m]); %I wanted to see what n=2 looks like

%sparse converts casadi.DM classtype to a list of doubles
plot(uVals,sparse(sVals)); 
hold on
%plot(uVals,sparse(zVals));
hold off

%generating normally distributed data points (~N(0,0.1))
disp(length(sparse(sVals)));
normRands = normrnd(0,0.1,[100,1]);
Redist_sVals = abs(normRands+sVals(:));%normally redistributed sVals
%scatter(uVals,sparse(Redist_sVals),'filled');

%generating more interesting normally distributed data points
%xNormRands = linspace(0,1,100); %preallocating array for efficiency
%x dependent random variables -> N(0,x*0.2)
%for i=1:100
    %xNormRands(i) = normrnd(0,0.2*sVals(i));
%end
%xRedist_sVals = abs(xNormRands+sVals(:));
%scatter(uVals,sparse(xRedist_sVals),'filled');

%properly generating normally distributed "experimental" data

inputVals_p = linspace(0,40,15);
outputVals_p = sparse(xStar(inputVals_p,[a,K,n]));



properVals_p = meshgrid(1:15,1:15);
for i=1:15
    properVals_p(:,i) = abs(outputVals_p(i)+normrnd(0,outputVals_p(i)*0.08,[15,1]));
end
disp(properVals_p);
hold on
for i=1:15
    scatter(ones(1,15)*inputVals_p(i),properVals_p(:,i));
end
hold off

%%
%%%stochastic simulation algorithm



