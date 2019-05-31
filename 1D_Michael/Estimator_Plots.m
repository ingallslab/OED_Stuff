Norm_D = dlmread('200x20__Normal_D_90.txt','\t');
Norm_Ds = dlmread('200x20__Normal_Ds_90.txt','\t');
Norm_L = dlmread('200x20__Normal_LinSpace_90.txt','\t');
SSA_D = dlmread('200x20__SSA_D_90.txt','\t');
SSA_Ds = dlmread('200x20__SSA_Ds_90.txt','\t');
SSA_L = dlmread('200x20__SSA_LinSpace_90.txt','\t');

% scatterEData(Norm_D,1);
% scatterEData(Norm_Ds,2);
% scatterEData(Norm_L,3);
% scatterEData(SSA_D,4);
% scatterEData(SSA_Ds,5);
% scatterEData(SSA_L,6);

h=figure(7);
set(h,'Position', [10 10 900 900]);
sz=25;
subplots=[];
ndm=[];
ndsm=[];
nlm=[];
sdm=[];
sdsm=[];
slm=[];

for i=1:4
    ndm=[ndm mean(Norm_D(:,i))];
    ndsm=[ndsm mean(Norm_Ds(:,i))];
    nlm=[nlm mean(Norm_L(:,i))];
    sdm=[sdm mean(SSA_D(:,i))];
    sdsm=[sdsm mean(SSA_Ds(:,i))];
    slm=[slm mean(SSA_L(:,i))];
end

theta_t = [0.5,3,9,3];
disp(norm(ndm-theta_t));
disp(norm(ndsm-theta_t));
disp(norm(nlm-theta_t));
disp(norm(sdm-theta_t));
disp(norm(sdsm-theta_t));
disp(norm(slm-theta_t));

subplots=[subplots subplot(3,3,1)];
hold on

plotErrorEllipse(subplots(1),ndm([1,2]),cov([Norm_D(:,1) Norm_D(:,2)]),[0, 255, 25]/255);
scatter(subplots(1),ndm(1),ndm(2),25,[0, 255, 25]/255,'filled');

plotErrorEllipse(subplots(1),ndsm([1,2]),cov([Norm_Ds(:,1) Norm_Ds(:,2)]),[0, 198, 19]/255);
scatter(subplots(1),ndsm(1),ndsm(2),25,[0, 198, 19]/255,'filled');

plotErrorEllipse(subplots(1),nlm([1,2]),cov([Norm_L(:,1) Norm_L(:,2)]),[0, 48, 0]/255);
scatter(subplots(1),nlm(1),nlm(2),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(1),sdm([1,2]),cov([SSA_D(:,1) SSA_D(:,2)]),[255, 0, 0]/255);
scatter(subplots(1),sdm(1),sdm(2),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(1),sdsm([1,2]),cov([SSA_Ds(:,1) SSA_Ds(:,2)]),[183, 0, 0]/255);
scatter(subplots(1),sdsm(1),sdsm(2),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(1),slm([1,2]),cov([SSA_L(:,1) SSA_L(:,2)]),[96, 0, 0]/255);
scatter(subplots(1),slm(1),slm(2),25,[96, 0, 0]/255,'filled');

scatter(subplots(1),0.5,3,40,'blue');

title(subplots(1),'\alpha versus \alpha_0');
xlabel(subplots(1),'\alpha_0'); ylabel(subplots(1),'\alpha');
hold off



subplots=[subplots subplot(3,3,2)];
hold on

plotErrorEllipse(subplots(2),ndm([1,3]),cov([Norm_D(:,1) Norm_D(:,3)]),[0, 255, 25]/255);
scatter(subplots(2),ndm(1),ndm(3),25,[0, 255, 25]/255,'filled');

plotErrorEllipse(subplots(2),ndsm([1,3]),cov([Norm_Ds(:,1) Norm_Ds(:,3)]),[0, 198, 19]/255);
scatter(subplots(2),ndsm(1),ndsm(3),25,[0, 198, 19]/255,'filled');

plotErrorEllipse(subplots(2),nlm([1,3]),cov([Norm_L(:,1) Norm_L(:,3)]),[0, 48, 0]/255);
scatter(subplots(2),nlm(1),nlm(3),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(2),sdm([1,3]),cov([SSA_D(:,1) SSA_D(:,3)]),[255, 0, 0]/255);
scatter(subplots(2),sdm(1),sdm(3),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(2),sdsm([1,3]),cov([SSA_Ds(:,1) SSA_Ds(:,3)]),[183, 0, 0]/255);
scatter(subplots(2),sdsm(1),sdsm(3),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(2),slm([1,3]),cov([SSA_L(:,1) SSA_L(:,3)]),[96, 0, 0]/255);
scatter(subplots(2),slm(1),slm(3),25,[96, 0, 0]/255,'filled');

scatter(subplots(2),0.5,9,40,'blue');

title(subplots(2),'K versus \alpha_0');
xlabel(subplots(2),'\alpha_0'); ylabel(subplots(2),'K');
hold off

subplots=[subplots subplot(3,3,3)];
hold on

plotErrorEllipse(subplots(3),ndm([1,4]),cov([Norm_D(:,1) Norm_D(:,4)]),[0, 255, 25]/255);
scatter(subplots(3),ndm(1),ndm(4),25,[0, 255, 25]/255,'filled');

plotErrorEllipse(subplots(3),ndsm([1,4]),cov([Norm_Ds(:,1) Norm_Ds(:,4)]),[0, 198, 19]/255);
scatter(subplots(3),ndsm(1),ndsm(4),25,[0, 198, 19]/255,'filled');

plotErrorEllipse(subplots(3),nlm([1,4]),cov([Norm_L(:,1) Norm_L(:,4)]),[0, 48, 0]/255);
scatter(subplots(3),nlm(1),nlm(4),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(3),sdm([1,4]),cov([SSA_D(:,1) SSA_D(:,4)]),[255, 0, 0]/255);
scatter(subplots(3),sdm(1),sdm(4),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(3),sdsm([1,4]),cov([SSA_Ds(:,1) SSA_Ds(:,4)]),[183, 0, 0]/255);
scatter(subplots(3),sdsm(1),sdsm(4),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(3),slm([1,4]),cov([SSA_L(:,1) SSA_L(:,4)]),[96, 0, 0]/255);
scatter(subplots(3),slm(1),slm(4),25,[96, 0, 0]/255,'filled');

scatter(subplots(3),0.5,3,40,'blue');

title(subplots(3),'n versus \alpha_0');
xlabel(subplots(3),'\alpha_0'); ylabel(subplots(3),'n');
hold off

subplots=[subplots subplot(3,3,4)];
hold on

plotErrorEllipse(subplots(4),ndm([2,3]),cov([Norm_D(:,2) Norm_D(:,3)]),[0, 255, 25]/255);
scatter(subplots(4),ndm(2),ndm(3),25,[0, 255, 25]/255,'filled');

plotErrorEllipse(subplots(4),ndsm([2,3]),cov([Norm_Ds(:,2) Norm_Ds(:,3)]),[0, 198, 19]/255);
scatter(subplots(4),ndsm(2),ndsm(3),25,[0, 198, 19]/255,'filled');

plotErrorEllipse(subplots(4),nlm([2,3]),cov([Norm_L(:,2) Norm_L(:,3)]),[0, 48, 0]/255);
scatter(subplots(4),nlm(2),nlm(3),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(4),sdm([2,3]),cov([SSA_D(:,2) SSA_D(:,3)]),[255, 0, 0]/255);
scatter(subplots(4),sdm(2),sdm(3),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(4),sdsm([2,3]),cov([SSA_Ds(:,2) SSA_Ds(:,3)]),[183, 0, 0]/255);
scatter(subplots(4),sdsm(2),sdsm(3),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(4),slm([2,3]),cov([SSA_L(:,2) SSA_L(:,3)]),[96, 0, 0]/255);
scatter(subplots(4),slm(2),slm(3),25,[96, 0, 0]/255,'filled');

scatter(subplots(4),3,9,40,'blue');

title(subplots(4),'K versus \alpha');
xlabel(subplots(4),'\alpha'); ylabel(subplots(4),'K');
hold off

subplots=[subplots subplot(3,3,5)];
hold on

plotErrorEllipse(subplots(5),ndm([2,4]),cov([Norm_D(:,2) Norm_D(:,4)]),[0, 255, 25]/255);
scatter(subplots(5),ndm(2),ndm(4),25,[0, 255, 25]/255,'filled');

plotErrorEllipse(subplots(5),ndsm([2,4]),cov([Norm_Ds(:,2) Norm_Ds(:,4)]),[0, 198, 19]/255);
scatter(subplots(5),ndsm(2),ndsm(4),25,[0, 198, 19]/255,'filled');

plotErrorEllipse(subplots(5),nlm([2,4]),cov([Norm_L(:,2) Norm_L(:,4)]),[0, 48, 0]/255);
scatter(subplots(5),nlm(2),nlm(4),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(5),sdm([2,4]),cov([SSA_D(:,2) SSA_D(:,4)]),[255, 0, 0]/255);
scatter(subplots(5),sdm(2),sdm(4),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(5),sdsm([2,4]),cov([SSA_Ds(:,2) SSA_Ds(:,4)]),[183, 0, 0]/255);
scatter(subplots(5),sdsm(2),sdsm(4),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(5),slm([2,4]),cov([SSA_L(:,2) SSA_L(:,4)]),[96, 0, 0]/255);
scatter(subplots(5),slm(2),slm(4),25,[96, 0, 0]/255,'filled');

scatter(subplots(5),3,3,40,'blue');

title(subplots(5),'n versus \alpha');
xlabel(subplots(5),'\alpha'); ylabel(subplots(5),'n');
hold off

subplots=[subplots subplot(3,3,6)];
hold on

plotErrorEllipse(subplots(6),ndm([3,4]),cov([Norm_D(:,3) Norm_D(:,4)]),[0, 255, 25]/255,'Normal D-Optimal');
scatter(subplots(6),ndm(3),ndm(4),25,[0, 255, 25]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),ndsm([3,4]),cov([Norm_Ds(:,3) Norm_Ds(:,4)]),[0, 198, 19]/255,'Normal Ds-Optimal');
scatter(subplots(6),ndsm(3),ndsm(4),25,[0, 198, 19]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),nlm([3,4]),cov([Norm_L(:,3) Norm_L(:,4)]),[0, 48, 0]/255,'Normal Non-OED');
scatter(subplots(6),nlm(3),nlm(4),25,[0, 48, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),sdm([3,4]),cov([SSA_D(:,3) SSA_D(:,4)]),[255, 0, 0]/255,'SSA D-Optimal');
scatter(subplots(6),sdm(3),sdm(4),25,[255, 0, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),sdsm([3,4]),cov([SSA_Ds(:,3) SSA_Ds(:,4)]),[183, 0, 0]/255,'SSA Ds-Optimal');
scatter(subplots(6),sdsm(3),sdsm(4),25,[183, 0, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),slm([3,4]),cov([SSA_L(:,3) SSA_L(:,4)]),[96, 0, 0]/255,'SSA Non-OED');
scatter(subplots(6),slm(3),slm(4),25,[96, 0, 0]/255,'filled','HandleVisibility','off');

scatter(subplots(6),9,3,40,'blue','HandleVisibility','off');

title(subplots(6),'n versus K');
xlabel(subplots(6),'K'); ylabel(subplots(6),'n');
hold off
hleg = subplot(3,3,7);
axis(hleg,'off');
pleg = get(hleg,'Position');
leg=legend(subplot(3,3,6));
set(leg,'Position',pleg);
saveas(h,'Covariance_Plot_200x20_90.png');

%%
% h = figure(3);
% h.Color = [1,1,1];
% 
% sz=25;
% 
% c1= subplot(2,3,1);
% hold on
% scatter(c1,ssaMeans(1),ssaMeans(2),40,'red','filled');
% plotErrorEllipse(c1,ssaMeans([1,2]),cov([SSAData(:,1) SSAData(:,2)]));
% 
% scatter(c1,normMeans(1),normMeans(2),40,'green','filled');
% plotErrorEllipse(c1,normMeans([1,2]),cov([NormalData(:,1) NormalData(:,2)]),'green');
% 
% scatter(c1,0.5,3,40,'blue','filled');
% hold off
% title(c1,'\alpha versus \alpha_0');
% xlabel(c1,'\alpha_0'); ylabel(c1,'\alpha');
% 
% c2=subplot(2,3,2);
% hold on
% scatter(c2,ssaMeans(1),ssaMeans(3),40,'red','filled');
% plotErrorEllipse(c2,ssaMeans([1,3]),cov([SSAData(:,1) SSAData(:,3)]));
% 
% scatter(c2,normMeans(1),normMeans(3),40,'green','filled');
% plotErrorEllipse(c2,normMeans([1,3]),cov([NormalData(:,1) NormalData(:,3)]),'green');
% 
% scatter(c2,0.5,9,40,'blue','filled');
% hold off
% title(c2,'K versus \alpha_0');
% xlabel(c2,'\alpha_0'); ylabel(c2,'K');
% 
% c3=subplot(2,3,3);
% hold on
% scatter(c3,ssaMeans(1),ssaMeans(4),40,'red','filled');
% plotErrorEllipse(c3,ssaMeans([1,4]),cov([SSAData(:,1) SSAData(:,4)]));
% 
% scatter(c3,normMeans(1),normMeans(4),40,'green','filled');
% plotErrorEllipse(c3,normMeans([1,4]),cov([NormalData(:,1) NormalData(:,4)]),'green');
% 
% scatter(c3,0.5,3,40,'blue','filled');
% hold off
% title(c3,'n versus \alpha_0');
% xlabel(c3,'\alpha_0'); ylabel(c3,'n');
% 
% c4=subplot(2,3,4);
% hold on
% scatter(c4,ssaMeans(2),ssaMeans(3),40,'red','filled');
% plotErrorEllipse(c4,ssaMeans([2,3]),cov([SSAData(:,2) SSAData(:,3)]));
% 
% scatter(c4,normMeans(2),normMeans(3),40,'green','filled');
% plotErrorEllipse(c4,normMeans([2,3]),cov([NormalData(:,2) NormalData(:,3)]),'green');
% 
% scatter(c4,3,9,40,'blue','filled');
% hold off
% title(c4,'K versus \alpha');
% xlabel(c4,'\alpha'); ylabel(c4,'K');
% 
% c5=subplot(2,3,5);
% hold on
% scatter(c5,ssaMeans(2),ssaMeans(4),40,'red','filled');
% plotErrorEllipse(c5,ssaMeans([2,4]),cov([SSAData(:,2) SSAData(:,4)]));
% 
% scatter(c5,normMeans(2),normMeans(4),40,'green','filled');
% plotErrorEllipse(c5,normMeans([2,4]),cov([NormalData(:,2) NormalData(:,4)]),'green');
% 
% scatter(c5,3,3,40,'blue','filled');
% hold off
% title(c5,'n versus \alpha');
% xlabel(c5,'\alpha'); ylabel(c5,'n');
% 
% c6=subplot(2,3,6);
% hold on
% scatter(c6,ssaMeans(3),ssaMeans(4),40,'red','filled');
% plotErrorEllipse(c6,ssaMeans([3,4]),cov([SSAData(:,3) SSAData(:,4)]));
% 
% scatter(c6,normMeans(3),normMeans(4),40,'green','filled');
% plotErrorEllipse(c6,normMeans([3,4]),cov([NormalData(:,3) NormalData(:,4)]),'green');
% 
function scatterEData(Data,fignum)
    means=[];
    for i=1:size(Data,2)
        means=[means mean(Data(:,i))];
    end
    figure(fignum);

    sz=25;

    b1= subplot(2,3,1);
    hold on
    scatter(b1,Data(:,1),Data(:,2),sz);
    scatter(b1,means(1),means(2),40,'red','filled');
    plotErrorEllipse(b1,means([1,2]),cov([Data(:,1) Data(:,2)]));
    hold off
    title(b1,'\alpha versus \alpha_0');
    xlabel(b1,'\alpha_0'); ylabel(b1,'\alpha');

    b2=subplot(2,3,2);
    hold on
    scatter(b2,Data(:,1),Data(:,3),sz);
    scatter(b2,means(1),means(3),40,'red','filled');
    plotErrorEllipse(b2,means([1,3]),cov([Data(:,1) Data(:,3)]));
    hold off
    title(b2,'K versus \alpha_0');
    xlabel(b2,'\alpha_0'); ylabel(b2,'K');

    b3=subplot(2,3,3);
    hold on
    scatter(b3,Data(:,1),Data(:,4),sz);
    scatter(b3,means(1),means(4),40,'red','filled');
    plotErrorEllipse(b3,means([1,4]),cov([Data(:,1) Data(:,4)]));
    hold off
    title(b3,'n versus \alpha_0');
    xlabel(b3,'\alpha_0'); ylabel(b3,'n');

    b4=subplot(2,3,4);
    hold on
    scatter(b4,Data(:,2),Data(:,3),sz);
    scatter(b4,means(2),means(3),40,'red','filled');
    plotErrorEllipse(b4,means([2,3]),cov([Data(:,2) Data(:,3)]));
    hold off
    title(b4,'K versus \alpha');
    xlabel(b4,'\alpha'); ylabel(b4,'K');

    b5=subplot(2,3,5);
    hold on
    scatter(b5,Data(:,2),Data(:,4),sz);
    scatter(b5,means(2),means(4),40,'red','filled');
    plotErrorEllipse(b5,means([2,4]),cov([Data(:,2) Data(:,4)]));
    hold off
    title(b5,'n versus \alpha');
    xlabel(b5,'\alpha'); ylabel(b5,'n');

    b6=subplot(2,3,6);
    hold on
    scatter(b6,Data(:,3),Data(:,4),sz);
    scatter(b6,means(3),means(4),40,'red','filled');
    plotErrorEllipse(b6,means([3,4]),cov([Data(:,3) Data(:,4)]));
    hold off
    title(b6,'K versus n');
    xlabel(b6,'n'); ylabel(b6,'K');
end

function plotErrorEllipse(fig,mu,Sigma,color,name)
    [V,D] = eig(Sigma);
    t=linspace(0,2*pi);
    a=(V*sqrt(D))*[cos(t(:))';sin(t(:))'];
    p=plot(fig,a(1,:)+mu(1),a(2,:)+mu(2));
    if exist('name','var')
        p.DisplayName=name;
    end
    p.LineWidth=1.5;
    if(exist('color','var'))
        p.Color=color;
    else
        p.Color='red';
    end
end
