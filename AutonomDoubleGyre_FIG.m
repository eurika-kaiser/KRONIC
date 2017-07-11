path2data = '../Data/'; 
load([path2data,[ModelName1,'Ensemble.mat']])
path2figs = '../Figures/DOUBLE_GYRE/'; mkdir(path2figs)

%% Phase plot with ensemble trajectories
cmap = gray(30); cmap=flipud(cmap(1:end,:));
figure,surf(x,y,StreamFun-10), shading interp, hold on, box on, view(2)
caxis([min(min(StreamFun))-10 max(max(StreamFun))-10])
colormap(cmap)

for iy = 1:Ny
    for ix = 1:Nx
        tmp = squeeze(yout_ctrl(:,iy,ix,:));
        [ynew,StartIDX] = applyPeriodicBC(tmp',[0,2],[0,1],Psi);
        for i = 1:length(StartIDX)-1
            plot(ynew(StartIDX(i):StartIDX(i+1)-1,1),ynew(StartIDX(i):StartIDX(i+1)-1,2),'-','Color',[.1,.1,1])
        end
        plot(ynew(1,1),ynew(1,2),'.','Color',[1 0 0],'MarkerSize',8)
    end
end
plot(ynew(1,1),ynew(1,2),'.','Color',[1 0 0],'MarkerSize',10)
contour(X,Y,StreamFun, [-1:0.05:1],'k','ShowText','on')
contour3(X,Y,StreamFun+20, [0.2+20, 0.2+20],'--','Color',0.8.*ones(1,3),... %[1,0,0]
    'ShowText','off','LineWidth',2)
% text(0.45,0.75,'0.2','FontSize',14,'Color',0.7.*ones(1,3))%[1,0,0]
box on
axis tight

xlabel('x1'), ylabel('x2')
set(gca,'FontSize',14, 'LineWidth',1)
axis equal
set(gcf,'Position',[100 100 600 350])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName1,'PhasePlot','.eps']);
