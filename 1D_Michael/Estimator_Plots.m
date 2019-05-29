NormalData = dlmread('./1D_Michael/Data/Norm_Linspace_LNA_90.txt','\t');
SSAData = dlmread('./1D_Michael/Data/SSA_Linspace_90.txt','\t');

normMeans=[];
for i=1:size(NormalData,2)
    normMeans=[normMeans mean(NormalData(:,i))];
end
ssaMeans=[];
for i=1:size(SSAData,2)
    ssaMeans=[ssaMeans mean(SSAData(:,i))];
end


g = figure(1);
g.Color = [1,1,1];

sz=25;

b1= subplot(2,3,1);
hold on
scatter(b1,NormalData(:,1),NormalData(:,2),sz);
scatter(b1,normMeans(1),normMeans(2),40,'red','filled');
plotErrorEllipse(b1,normMeans([1,2]),cov([NormalData(:,1) NormalData(:,2)]));
hold off
title(b1,'\alpha versus \alpha_0');
xlabel(b1,'\alpha_0'); ylabel(b1,'\alpha');

b2=subplot(2,3,2);
hold on
scatter(b2,NormalData(:,1),NormalData(:,3),sz);
scatter(b2,normMeans(1),normMeans(3),40,'red','filled');
plotErrorEllipse(b2,normMeans([1,3]),cov([NormalData(:,1) NormalData(:,3)]));
hold off
title(b2,'K versus \alpha_0');
xlabel(b2,'\alpha_0'); ylabel(b2,'K');

b3=subplot(2,3,3);
hold on
scatter(b3,NormalData(:,1),NormalData(:,4),sz);
scatter(b3,normMeans(1),normMeans(4),40,'red','filled');
plotErrorEllipse(b3,normMeans([1,4]),cov([NormalData(:,1) NormalData(:,4)]));
hold off
title(b3,'n versus \alpha_0');
xlabel(b3,'\alpha_0'); ylabel(b3,'n');

b4=subplot(2,3,4);
hold on
scatter(b4,NormalData(:,2),NormalData(:,3),sz);
scatter(b4,normMeans(2),normMeans(3),40,'red','filled');
plotErrorEllipse(b4,normMeans([2,3]),cov([NormalData(:,2) NormalData(:,3)]));
hold off
title(b4,'K versus \alpha');
xlabel(b4,'\alpha'); ylabel(b4,'K');

b5=subplot(2,3,5);
hold on
scatter(b5,NormalData(:,2),NormalData(:,4),sz);
scatter(b5,normMeans(2),normMeans(4),40,'red','filled');
plotErrorEllipse(b5,normMeans([2,4]),cov([NormalData(:,2) NormalData(:,4)]));
hold off
title(b5,'n versus \alpha');
xlabel(b5,'\alpha'); ylabel(b5,'n');

b6=subplot(2,3,6);
hold on
scatter(b6,NormalData(:,3),NormalData(:,4),sz);
scatter(b6,normMeans(3),normMeans(4),40,'red','filled');
plotErrorEllipse(b6,normMeans([3,4]),cov([NormalData(:,3) NormalData(:,4)]));
hold off
title(b6,'K versus n');
xlabel(b6,'n'); ylabel(b6,'K');

f = figure(2);
f.Color = [1,1,1];

sz=25;

a1= subplot(2,3,1);
hold on
scatter(a1,SSAData(:,1),SSAData(:,2),sz);
scatter(a1,ssaMeans(1),ssaMeans(2),40,'red','filled');
plotErrorEllipse(a1,ssaMeans([1,2]),cov([SSAData(:,1) SSAData(:,2)]));
hold off
title(a1,'\alpha versus \alpha_0');
xlabel(a1,'\alpha_0'); ylabel(a1,'\alpha');

a2=subplot(2,3,2);
hold on
scatter(a2,SSAData(:,1),SSAData(:,3),sz);
scatter(a2,ssaMeans(1),ssaMeans(3),40,'red','filled');
plotErrorEllipse(a2,ssaMeans([1,3]),cov([SSAData(:,1) SSAData(:,3)]));
hold off
title(a2,'K versus \alpha_0');
xlabel(a2,'\alpha_0'); ylabel(a2,'K');

a3=subplot(2,3,3);
hold on
scatter(a3,SSAData(:,1),SSAData(:,4),sz);
scatter(a3,ssaMeans(1),ssaMeans(4),40,'red','filled');
plotErrorEllipse(a3,ssaMeans([1,4]),cov([SSAData(:,1) SSAData(:,4)]));
hold off
title(a3,'n versus \alpha_0');
xlabel(a3,'\alpha_0'); ylabel(a3,'n');

a4=subplot(2,3,4);
hold on
scatter(a4,SSAData(:,2),SSAData(:,3),sz);
scatter(a4,ssaMeans(2),ssaMeans(3),40,'red','filled');
plotErrorEllipse(a4,ssaMeans([2,3]),cov([SSAData(:,2) SSAData(:,3)]));
hold off
title(a4,'K versus \alpha');
xlabel(a4,'\alpha'); ylabel(a4,'K');

a5=subplot(2,3,5);
hold on
scatter(a5,SSAData(:,2),SSAData(:,4),sz);
scatter(a5,ssaMeans(2),ssaMeans(4),40,'red','filled');
plotErrorEllipse(a5,ssaMeans([2,4]),cov([SSAData(:,2) SSAData(:,4)]));
hold off
title(a5,'n versus \alpha');
xlabel(a5,'\alpha'); ylabel(a5,'n');

a6=subplot(2,3,6);
hold on
scatter(a6,SSAData(:,3),SSAData(:,4),sz);
scatter(a6,ssaMeans(3),ssaMeans(4),40,'red','filled');
plotErrorEllipse(a6,ssaMeans([3,4]),cov([SSAData(:,3) SSAData(:,4)]));
hold off
title(a6,'K versus n');
xlabel(a6,'n'); ylabel(a6,'K');

%%
h = figure(3);
h.Color = [1,1,1];

sz=25;

c1= subplot(2,3,1);
hold on
scatter(c1,ssaMeans(1),ssaMeans(2),40,'red','filled');
plotErrorEllipse(c1,ssaMeans([1,2]),cov([SSAData(:,1) SSAData(:,2)]));

scatter(c1,normMeans(1),normMeans(2),40,'green','filled');
plotErrorEllipse(c1,normMeans([1,2]),cov([NormalData(:,1) NormalData(:,2)]),'green');

scatter(c1,0.5,3,40,'blue','filled');
hold off
title(c1,'\alpha versus \alpha_0');
xlabel(c1,'\alpha_0'); ylabel(c1,'\alpha');

c2=subplot(2,3,2);
hold on
scatter(c2,ssaMeans(1),ssaMeans(3),40,'red','filled');
plotErrorEllipse(c2,ssaMeans([1,3]),cov([SSAData(:,1) SSAData(:,3)]));

scatter(c2,normMeans(1),normMeans(3),40,'green','filled');
plotErrorEllipse(c2,normMeans([1,3]),cov([NormalData(:,1) NormalData(:,3)]),'green');

scatter(c2,0.5,9,40,'blue','filled');
hold off
title(c2,'K versus \alpha_0');
xlabel(c2,'\alpha_0'); ylabel(c2,'K');

c3=subplot(2,3,3);
hold on
scatter(c3,ssaMeans(1),ssaMeans(4),40,'red','filled');
plotErrorEllipse(c3,ssaMeans([1,4]),cov([SSAData(:,1) SSAData(:,4)]));

scatter(c3,normMeans(1),normMeans(4),40,'green','filled');
plotErrorEllipse(c3,normMeans([1,4]),cov([NormalData(:,1) NormalData(:,4)]),'green');

scatter(c3,0.5,3,40,'blue','filled');
hold off
title(c3,'n versus \alpha_0');
xlabel(c3,'\alpha_0'); ylabel(c3,'n');

c4=subplot(2,3,4);
hold on
scatter(c4,ssaMeans(2),ssaMeans(3),40,'red','filled');
plotErrorEllipse(c4,ssaMeans([2,3]),cov([SSAData(:,2) SSAData(:,3)]));

scatter(c4,normMeans(2),normMeans(3),40,'green','filled');
plotErrorEllipse(c4,normMeans([2,3]),cov([NormalData(:,2) NormalData(:,3)]),'green');

scatter(c4,3,9,40,'blue','filled');
hold off
title(c4,'K versus \alpha');
xlabel(c4,'\alpha'); ylabel(c4,'K');

c5=subplot(2,3,5);
hold on
scatter(c5,ssaMeans(2),ssaMeans(4),40,'red','filled');
plotErrorEllipse(c5,ssaMeans([2,4]),cov([SSAData(:,2) SSAData(:,4)]));

scatter(c5,normMeans(2),normMeans(4),40,'green','filled');
plotErrorEllipse(c5,normMeans([2,4]),cov([NormalData(:,2) NormalData(:,4)]),'green');

scatter(c5,3,3,40,'blue','filled');
hold off
title(c5,'n versus \alpha');
xlabel(c5,'\alpha'); ylabel(c5,'n');

c6=subplot(2,3,6);
hold on
scatter(c6,ssaMeans(3),ssaMeans(4),40,'red','filled');
plotErrorEllipse(c6,ssaMeans([3,4]),cov([SSAData(:,3) SSAData(:,4)]));

scatter(c6,normMeans(3),normMeans(4),40,'green','filled');
plotErrorEllipse(c6,normMeans([3,4]),cov([NormalData(:,3) NormalData(:,4)]),'green');

scatter(c6,9,3,40,'blue','filled');
hold off
title(c6,'K versus n');
xlabel(c6,'n'); ylabel(c6,'K');




function plotErrorEllipse(fig,mu,Sigma,color)
    disp(Sigma);
    [V,D] = eig(Sigma);
    t=linspace(0,2*pi);
    a=(V*sqrt(D))*[cos(t(:))';sin(t(:))'];
    p=plot(fig,a(1,:)+mu(1),a(2,:)+mu(2));
    p.LineWidth=1.5;
    if(exist('color','var'))
        p.Color=color;
    else
        p.Color='red';
    end
end
