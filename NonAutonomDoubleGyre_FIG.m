
path2data = '../Data/'; 
load([path2data,[ModelName1,'Data.mat']])
path2figs = '../Figures/DOUBLE_GYRE/'; mkdir(path2figs)

%% Phase plot with single trajectory
clear ph
[y1_corrected,StartIDX,outside_vec,switched_domain_vec,Hvec] = applyPeriodicBC(y1,[0,2],[0,1],Hsteady);
cmap = gray(30); cmap=flipud(cmap(1:end,:));

figure,surf(x,y,StreamFun-10), shading interp, hold on, box on, view(2)
caxis([min(min(StreamFun))-10 max(max(StreamFun))-10])
colormap(cmap)
caxis([min(min(StreamFun))-10 max(max(StreamFun))-10])
contour(X,Y,StreamFun, [-1:0.05:1],'k','ShowText','on')

ph(1) = plot(y0(:,1),y0(:,2),'-','Color',[0.8,0.8,0],'LineWidth',2);
plot(y1_corrected(1,1),y1_corrected(1,2),'.','Color',[1 0 0],'MarkerSize',20)
for i = 1:length(StartIDX)-1
    ph(2) = plot(y1_corrected(StartIDX(i):StartIDX(i+1)-1,1),y1_corrected(StartIDX(i):StartIDX(i+1)-1,2),'-','Color',[.1,.1,1],'LineWidth',2);
end
contour3(X,Y,PSIfield(:,:,end)+20, [0.2+20, 0.2+20],'--','Color',0.8.*ones(1,3), ...
    'ShowText','off','LineWidth',2)
% text(0.45,0.75,'-0.2','FontSize',14,'Color',[0.8,0.8,1])
box on
% legend(ph,'Unforced', 'Controlled')
xlabel('x1'), ylabel('x2')
set(gca,'FontSize',14, 'LineWidth',1)
axis equal
set(gcf,'Position',[100 100 600 350])
set(gcf,'PaperPositionMode','auto')
print('-opengl','-depsc2', '-loose', [path2figs,ModelName1,'singleIC_PhasePlot','.eps']);

%% Cost etc. for single trajectory
LineWidth = 2;
fhandle = figure; hold on, box on, grid on
plot(tspan,Psivals0,'-','Color',[0.8,0.8,0],'LineWidth',2);
plot(tspan,Psivals1,'-','Color',[.1,.1,1],'LineWidth',2)
xlabel('t'), ylabel('Psi')
set(gca,'FontSize',16)
axis tight
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'singleIC_Psi','.eps']);

clear ph
fhandle = figure; hold on, box on, grid on
ph(1) = plot(tspan,Jvals0,'-','Color',[0.8,0.8,0],'LineWidth',2);
ph(2) = plot(tspan,Jvals1,'-','Color',[.1,.1,1],'LineWidth',2);
xlabel('t'), ylabel('J')
set(gca,'FontSize',16)
axis tight
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'singleIC_J','.eps']);

figure; hold on
ph(1) = plot([0,1],[0,0],'-','Color',[0.8,0.8,0],'LineWidth',6);
ph(2) = plot([0,1],[0,0],'-','Color',[.1,.1,1],'LineWidth',6);
pl = legend('Unforced', 'Controlled'); pl.Position = [0.3 0.35 0.5 0.2];
xlim([-1 2]), axis off
set(gca,'FontSize',16)
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'singleIC_legend','.eps']);


clear ph
fhandle = figure; hold on, box on, grid on
% plot(tspan,uvals0,'-','Color',[0.8,0.8,0],'LineWidth',2);
ph(1) = plot(tspan,uvals1(:,1),'-','Color',[.1,.1,1],'LineWidth',2);
ph(2) = plot(tspan,uvals1(:,2),'-','Color',[.1,0.7,.1],'LineWidth',2);
xlabel('t'), ylabel('u')
set(gca,'FontSize',16)
axis tight
legend('u1','u2')
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'singleIC_u','.eps']);
