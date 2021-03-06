Norm_D = dlmread('500x8__Normal_D_90.txt','\t');
Norm_Ds = dlmread('500x8__Normal_Ds_90.txt','\t');
Norm_L = dlmread('500x8__Normal_LinSpace_90.txt','\t');
SSA_D = dlmread('500x8__SSA_D_90.txt','\t');
SSA_Ds = dlmread('500x8__SSA_Ds_90.txt','\t');
SSA_L = dlmread('500x8__SSA_LinSpace_90.txt','\t');

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


Omega = 90;
numExp=20;
u = [ linspace(0,0.1,32) linspace(0.1,0.2,64) linspace(0.2,0.3,32)];
inc = length(u)/numExp;
if inc<1
    disp('Too many experiments!');
    return
end 
u = u(1:round(inc):end);  



co=inv(sum(200*FIM_comp([0.5,3,9,3,-18,117],Omega,u,[0.1,0.2]),3));
co=co(1:4,1:4);
co12 = co(1:2,1:2);
co13 = [[co(1,1) co(1,3)];[co(3,1),co(3,3)]];
co14 = [[co(1,1) co(1,4)];[co(4,1),co(4,4)]];
co23 = [[co(2,2) co(2,3)];[co(3,2),co(3,3)]];
co24 = [[co(2,2) co(2,4)];[co(4,2),co(4,4)]];
co34 = [[co(3,3) co(3,4)];[co(4,3),co(4,4)]];

[u_opt, w_opt]=D_opt_c([0.5,3,9,3,-18,117],Omega,u,FIM_comp([0.5,3,9,3,-18,117],Omega,u,[0.1,0.2]),[0.1,0.2]);
numSamp=round(w_opt*200*20);
fc = FIM_comp([0.5,3,9,3,-18,117],Omega,u_opt,[0.1,0.2]);
Co=0;
for i=1:length(u_opt)
    Co = Co + numSamp(i)*fc(:,:,i);
end
Co=Co(1:4,1:4);
Co=inv(Co);

Co12 = Co(1:2,1:2);
Co13 = [[Co(1,1) Co(1,3)];[Co(3,1),Co(3,3)]];
Co14 = [[Co(1,1) Co(1,4)];[Co(4,1),Co(4,4)]];
Co23 = [[Co(2,2) Co(2,3)];[Co(3,2),Co(3,3)]];
Co24 = [[Co(2,2) Co(2,4)];[Co(4,2),Co(4,4)]];
Co34 = [[Co(3,3) Co(3,4)];[Co(4,3),Co(4,4)]];

[u2_opt, w2_opt]=Ds_opt([0.5,3,9,3,-18,117],Omega,u,FIM_comp([0.5,3,9,3,-18,117],Omega,u,[0.1,0.2]),[0.1,0.2]);
numSamp=round(w2_opt*200*20);
fc = FIM_comp([0.5,3,9,3,-18,117],Omega,u2_opt,[0.1,0.2]);
Co2=0;
for i=1:length(u2_opt)
    Co2 = Co2 + numSamp(i)*fc(:,:,i);
end
Co2=Co2(1:4,1:4);
Co2=inv(Co2);
CO12 = Co2(1:2,1:2);
CO13 = [[Co2(1,1) Co2(1,3)];[Co2(3,1),Co2(3,3)]];
CO14 = [[Co2(1,1) Co2(1,4)];[Co2(4,1),Co2(4,4)]];
CO23 = [[Co2(2,2) Co2(2,3)];[Co2(3,2),Co2(3,3)]];
CO24 = [[Co2(2,2) Co2(2,4)];[Co2(4,2),Co2(4,4)]];
CO34 = [[Co2(3,3) Co2(3,4)];[Co2(4,3),Co2(4,4)]];



subplots=[subplots subplot(3,3,1)];
hold on
% 
% plotErrorEllipse(subplots(1),ndm([1,2]),cov([Norm_D(:,1) Norm_D(:,2)]),[0, 255, 25]/255);
% scatter(subplots(1),ndm(1),ndm(2),25,[0, 255, 25]/255,'filled');
% 
% plotErrorEllipse(subplots(1),ndsm([1,2]),cov([Norm_Ds(:,1) Norm_Ds(:,2)]),[0, 198, 19]/255);
% scatter(subplots(1),ndsm(1),ndsm(2),25,[0, 198, 19]/255,'filled');
% 
% plotErrorEllipse(subplots(1),nlm([1,2]),cov([Norm_L(:,1) Norm_L(:,2)]),[0, 48, 0]/255);
% scatter(subplots(1),nlm(1),nlm(2),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(1),sdm([1,2]),cov([SSA_D(:,1) SSA_D(:,2)]),[255, 0, 0]/255);
scatter(subplots(1),sdm(1),sdm(2),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(1),sdsm([1,2]),cov([SSA_Ds(:,1) SSA_Ds(:,2)]),[183, 0, 0]/255);
scatter(subplots(1),sdsm(1),sdsm(2),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(1),slm([1,2]),cov([SSA_L(:,1) SSA_L(:,2)]),[96, 0, 0]/255);
scatter(subplots(1),slm(1),slm(2),25,[96, 0, 0]/255,'filled');

plotErrorEllipse(subplots(1),[0.5,3],CO12,[0, 20, 53]/255);
plotErrorEllipse(subplots(1),[0.5,3],co12,[2, 56, 142]/255);
plotErrorEllipse(subplots(1),[0.5,3],Co12,'blue');
scatter(subplots(1),0.5,3,40,'blue');

title(subplots(1),'\alpha versus \alpha_0');
xlabel(subplots(1),'\alpha_0'); ylabel(subplots(1),'\alpha');
hold off

subplots=[subplots subplot(3,3,2)];
hold on

% plotErrorEllipse(subplots(2),ndm([1,3]),cov([Norm_D(:,1) Norm_D(:,3)]),[0, 255, 25]/255);
% scatter(subplots(2),ndm(1),ndm(3),25,[0, 255, 25]/255,'filled');
% 
% plotErrorEllipse(subplots(2),ndsm([1,3]),cov([Norm_Ds(:,1) Norm_Ds(:,3)]),[0, 198, 19]/255);
% scatter(subplots(2),ndsm(1),ndsm(3),25,[0, 198, 19]/255,'filled');
% 
% plotErrorEllipse(subplots(2),nlm([1,3]),cov([Norm_L(:,1) Norm_L(:,3)]),[0, 48, 0]/255);
% scatter(subplots(2),nlm(1),nlm(3),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(2),sdm([1,3]),cov([SSA_D(:,1) SSA_D(:,3)]),[255, 0, 0]/255);
scatter(subplots(2),sdm(1),sdm(3),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(2),sdsm([1,3]),cov([SSA_Ds(:,1) SSA_Ds(:,3)]),[183, 0, 0]/255);
scatter(subplots(2),sdsm(1),sdsm(3),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(2),slm([1,3]),cov([SSA_L(:,1) SSA_L(:,3)]),[96, 0, 0]/255);
scatter(subplots(2),slm(1),slm(3),25,[96, 0, 0]/255,'filled');

plotErrorEllipse(subplots(2),[0.5,9],CO13,[0, 20, 53]/255);
plotErrorEllipse(subplots(2),[0.5,9],co13,[2, 56, 142]/255);
plotErrorEllipse(subplots(2),[0.5,9],Co13,'blue');
scatter(subplots(2),0.5,9,40,'blue');

title(subplots(2),'K versus \alpha_0');
xlabel(subplots(2),'\alpha_0'); ylabel(subplots(2),'K');
hold off

subplots=[subplots subplot(3,3,3)];
hold on
% 
% plotErrorEllipse(subplots(3),ndm([1,4]),cov([Norm_D(:,1) Norm_D(:,4)]),[0, 255, 25]/255);
% scatter(subplots(3),ndm(1),ndm(4),25,[0, 255, 25]/255,'filled');
% 
% plotErrorEllipse(subplots(3),ndsm([1,4]),cov([Norm_Ds(:,1) Norm_Ds(:,4)]),[0, 198, 19]/255);
% scatter(subplots(3),ndsm(1),ndsm(4),25,[0, 198, 19]/255,'filled');
% 
% plotErrorEllipse(subplots(3),nlm([1,4]),cov([Norm_L(:,1) Norm_L(:,4)]),[0, 48, 0]/255);
% scatter(subplots(3),nlm(1),nlm(4),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(3),sdm([1,4]),cov([SSA_D(:,1) SSA_D(:,4)]),[255, 0, 0]/255);
scatter(subplots(3),sdm(1),sdm(4),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(3),sdsm([1,4]),cov([SSA_Ds(:,1) SSA_Ds(:,4)]),[183, 0, 0]/255);
scatter(subplots(3),sdsm(1),sdsm(4),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(3),slm([1,4]),cov([SSA_L(:,1) SSA_L(:,4)]),[96, 0, 0]/255);
scatter(subplots(3),slm(1),slm(4),25,[96, 0, 0]/255,'filled');

plotErrorEllipse(subplots(3),[0.5,3],CO14,[0, 20, 53]/255);
plotErrorEllipse(subplots(3),[0.5,3],co14,[2, 56, 142]/255);
plotErrorEllipse(subplots(3),[0.5,3],Co14,'blue');
scatter(subplots(3),0.5,3,40,'blue');

title(subplots(3),'n versus \alpha_0');
xlabel(subplots(3),'\alpha_0'); ylabel(subplots(3),'n');
hold off

subplots=[subplots subplot(3,3,4)];
hold on

% plotErrorEllipse(subplots(4),ndm([2,3]),cov([Norm_D(:,2) Norm_D(:,3)]),[0, 255, 25]/255);
% scatter(subplots(4),ndm(2),ndm(3),25,[0, 255, 25]/255,'filled');
% 
% plotErrorEllipse(subplots(4),ndsm([2,3]),cov([Norm_Ds(:,2) Norm_Ds(:,3)]),[0, 198, 19]/255);
% scatter(subplots(4),ndsm(2),ndsm(3),25,[0, 198, 19]/255,'filled');
% 
% plotErrorEllipse(subplots(4),nlm([2,3]),cov([Norm_L(:,2) Norm_L(:,3)]),[0, 48, 0]/255);
% scatter(subplots(4),nlm(2),nlm(3),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(4),sdm([2,3]),cov([SSA_D(:,2) SSA_D(:,3)]),[255, 0, 0]/255);
scatter(subplots(4),sdm(2),sdm(3),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(4),sdsm([2,3]),cov([SSA_Ds(:,2) SSA_Ds(:,3)]),[183, 0, 0]/255);
scatter(subplots(4),sdsm(2),sdsm(3),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(4),slm([2,3]),cov([SSA_L(:,2) SSA_L(:,3)]),[96, 0, 0]/255);
scatter(subplots(4),slm(2),slm(3),25,[96, 0, 0]/255,'filled');

plotErrorEllipse(subplots(4),[3,9],CO23,[0, 20, 53]/255);
plotErrorEllipse(subplots(4),[3,9],co23,[2, 56, 142]/255);
plotErrorEllipse(subplots(4),[3,9],Co23,'blue');
scatter(subplots(4),3,9,40,'blue');

title(subplots(4),'K versus \alpha');
xlabel(subplots(4),'\alpha'); ylabel(subplots(4),'K');
hold off

subplots=[subplots subplot(3,3,5)];
hold on
% 
% plotErrorEllipse(subplots(5),ndm([2,4]),cov([Norm_D(:,2) Norm_D(:,4)]),[0, 255, 25]/255);
% scatter(subplots(5),ndm(2),ndm(4),25,[0, 255, 25]/255,'filled');
% 
% plotErrorEllipse(subplots(5),ndsm([2,4]),cov([Norm_Ds(:,2) Norm_Ds(:,4)]),[0, 198, 19]/255);
% scatter(subplots(5),ndsm(2),ndsm(4),25,[0, 198, 19]/255,'filled');
% 
% plotErrorEllipse(subplots(5),nlm([2,4]),cov([Norm_L(:,2) Norm_L(:,4)]),[0, 48, 0]/255);
% scatter(subplots(5),nlm(2),nlm(4),25,[0, 48, 0]/255,'filled');

plotErrorEllipse(subplots(5),sdm([2,4]),cov([SSA_D(:,2) SSA_D(:,4)]),[255, 0, 0]/255);
scatter(subplots(5),sdm(2),sdm(4),25,[255, 0, 0]/255,'filled');

plotErrorEllipse(subplots(5),sdsm([2,4]),cov([SSA_Ds(:,2) SSA_Ds(:,4)]),[183, 0, 0]/255);
scatter(subplots(5),sdsm(2),sdsm(4),25,[183, 0, 0]/255,'filled');

plotErrorEllipse(subplots(5),slm([2,4]),cov([SSA_L(:,2) SSA_L(:,4)]),[96, 0, 0]/255);
scatter(subplots(5),slm(2),slm(4),25,[96, 0, 0]/255,'filled');

plotErrorEllipse(subplots(5),[3,3],CO24,[0, 20, 53]/255);
plotErrorEllipse(subplots(5),[3,3],co24,[2, 56, 142]/255);
plotErrorEllipse(subplots(5),[3,3],Co24,'blue');
scatter(subplots(5),3,3,40,'blue');

title(subplots(5),'n versus \alpha');
xlabel(subplots(5),'\alpha'); ylabel(subplots(5),'n');
hold off

subplots=[subplots subplot(3,3,6)];
hold on
% 
% plotErrorEllipse(subplots(6),ndm([3,4]),cov([Norm_D(:,3) Norm_D(:,4)]),[0, 255, 25]/255,'Normal D-Optimal');
% scatter(subplots(6),ndm(3),ndm(4),25,[0, 255, 25]/255,'filled','HandleVisibility','off');
% 
% plotErrorEllipse(subplots(6),ndsm([3,4]),cov([Norm_Ds(:,3) Norm_Ds(:,4)]),[0, 198, 19]/255,'Normal Ds-Optimal');
% scatter(subplots(6),ndsm(3),ndsm(4),25,[0, 198, 19]/255,'filled','HandleVisibility','off');
% 
% plotErrorEllipse(subplots(6),nlm([3,4]),cov([Norm_L(:,3) Norm_L(:,4)]),[0, 48, 0]/255,'Normal Non-OED');
% scatter(subplots(6),nlm(3),nlm(4),25,[0, 48, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),sdm([3,4]),cov([SSA_D(:,3) SSA_D(:,4)]),[255, 0, 0]/255,'SSA D-Optimal');
scatter(subplots(6),sdm(3),sdm(4),25,[255, 0, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),sdsm([3,4]),cov([SSA_Ds(:,3) SSA_Ds(:,4)]),[183, 0, 0]/255,'SSA Ds-Optimal');
scatter(subplots(6),sdsm(3),sdsm(4),25,[183, 0, 0]/255,'filled','HandleVisibility','off');

plotErrorEllipse(subplots(6),slm([3,4]),cov([SSA_L(:,3) SSA_L(:,4)]),[96, 0, 0]/255,'SSA Non-OED');
scatter(subplots(6),slm(3),slm(4),25,[96, 0, 0]/255,'filled','HandleVisibility','off');


plotErrorEllipse(subplots(6),[9,3],CO34,[0, 20, 53]/255, 'FIM Ds-Optimal');
plotErrorEllipse(subplots(6),[9,3],co34,[2, 56, 142]/255, 'FIM Suboptimal');
plotErrorEllipse(subplots(6),[9,3],Co34,'blue','FIM D-Optimal');
scatter(subplots(6),9,3,40,'blue','HandleVisibility','off');

title(subplots(6),'n versus K');
xlabel(subplots(6),'K'); ylabel(subplots(6),'n');
hold off

leg=legend(subplots(6));
set(leg,'Location','southeast')
saveas(h,'Covariance_Plot_500x8_90.png');

% figure
% [X,Y]=meshgrid(8.4:0.001:9.2,2.65:0.001:3.05);
% qform1=zeros(size(X,1),size(X,2));
% qform2=zeros(size(X,1),size(X,2));
% qform3=zeros(size(X,1),size(X,2));
% for i=1:size(X,1)
%     for j=1:size(Y,2)
%         w = [9-X(i,j);3-Y(i,j)];
%         qform1(i,j) = w'*(co34)*w;
%         qform2(i,j) = w'*(Co34)*w;
%         qform3(i,j) = w'*(CO34)*w;
%     end
% end
% hold on
% s=surf(X,Y,qform1);
% s.EdgeColor = 'none';
% %surf(X,Y,qform2,[1,1],'blue');
% %surf(X,Y,qform3,[1,1],'blue');
% hold off
disp('Biases');
disp(norm(ndm-theta_t));
disp(norm(ndsm-theta_t));
disp(norm(nlm-theta_t));
disp(norm(sdm-theta_t));
disp(norm(sdsm-theta_t));
disp(norm(slm-theta_t));
disp('Determinants');
disp(det(cov(Norm_D)));
disp(det(cov(Norm_Ds)));
disp(det(cov(Norm_L)));
disp(det(cov(SSA_D)));
disp(det(cov(SSA_Ds)));
disp(det(cov(SSA_L)));
disp('Information');
disp(det(co));
disp(det(Co));
disp(det(Co2));



disp('MSE');
sMSEd = [0 0 0 0];
for i=1:size(SSA_D,1)
    sMSEd=sMSEd+(theta_t - SSA_D(i,:)).*(theta_t - SSA_D(i,:))/size(SSA_D,1);
end
disp(norm(sMSEd));
sMSEds = [0 0 0 0];
for i=1:size(SSA_D,1)
    sMSEds=sMSEds+(theta_t - SSA_Ds(i,:)).*(theta_t - SSA_Ds(i,:))/size(SSA_Ds,1);
end
disp(norm(sMSEds));
sMSEl = [0 0 0 0];
for i=1:size(SSA_L,1)
    sMSEl=sMSEl+(theta_t - SSA_L(i,:)).*(theta_t - SSA_L(i,:))/size(SSA_L,1);
end
disp(norm(sMSEl));

nMSEd = [0 0 0 0];
for i=1:size(Norm_D,1)
    nMSEd=nMSEd+(theta_t - Norm_D(i,:)).*(theta_t - Norm_D(i,:))/size(Norm_D,1);
end
disp(norm(nMSEd));
nMSEds = [0 0 0 0];
for i=1:size(Norm_D,1)
    nMSEds=nMSEds+(theta_t - Norm_Ds(i,:)).*(theta_t - Norm_Ds(i,:))/size(Norm_Ds,1);
end
disp(norm(nMSEds));
nMSEl = [0 0 0 0];
for i=1:size(Norm_L,1)
    nMSEl=nMSEl+(theta_t - Norm_L(i,:)).*(theta_t - Norm_L(i,:))/size(Norm_L,1);
end
disp(norm(nMSEl));

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
