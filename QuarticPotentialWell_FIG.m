path2data = '../Data/'; 
load([path2data,[ModelName1,'Data.mat']])
path2figs = '../Figures/QUARTIC_POT_WELL/'; mkdir(path2figs)

%% Show phase plot
cmap = jet(10);
cmap = cmap(1:2:end,:);

figure;
x = -0.88:0.01:0.88; hold on, box on
% y = +/- sqrt(2*REF- 1/2 *x^4)
trihandle = trisurf(tri, XY(:,1), XY(:,2),(Hxy)); hold on, shading interp, 
caxis([0 2])
colormap(gray(8))
view(2)
plot3(y1(:,1),y1(:,2),(max(Hxy)+1).*ones(size(y1(:,1))),'-r', 'LineWidth',2,'Color',[0,0.5,1])
plot3(y2(:,1),y2(:,2),(max(Hxy)+1).*ones(size(y2(:,1))),'-', 'LineWidth',2,'Color',[0,0.5,1])
plot3([-1.6 1.6],[-1.6 -1.6],[(max(Hxy)+1) (max(Hxy)+1)],'-k', 'LineWidth',1) %b
plot3([-1.6 1.6],[1.6 1.6],[(max(Hxy)+1) (max(Hxy)+1)],'-k', 'LineWidth',1) %t
plot3([-1.6 -1.6],[-1.6 1.6],[(max(Hxy)+1) (max(Hxy)+1)],'-k', 'LineWidth',1) %l
plot3([1.6 1.6],[-1.6 1.6],[(max(Hxy)+1) (max(Hxy)+1)],'-k', 'LineWidth',1) %r
plot3(x,sqrt(2*REF- 1/2 *x.^4),(max(Hxy)+1).*ones(size(x)),'--', 'LineWidth',2,'Color',[0.9 0.9 0.9])
plot3(x,-sqrt(2*REF- 1/2 *x.^4),(max(Hxy)+1).*ones(size(x)),'--', 'LineWidth',2,'Color',[0.9 0.9 0.9])
xlabel('x_1'), ylabel('x_2')
axis equal
axis([-1.6 1.6 -1.6 1.6])
set(gca,'FontSize',16,'LineWidth',1,'ytick',[-1 0 1])
set(gcf,'Position',[100 100 300 300])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [figpath,ModelName1,'PhasePlot','.eps']);