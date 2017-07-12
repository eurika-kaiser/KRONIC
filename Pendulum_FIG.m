path2data = '../Data/'; ModelName1='Pendulum_B10_';
load([path2data,[ModelName1,'Data.mat']])
path2figs = '../Figures/PENDULUM/'; mkdir(path2figs)


%% Show phase plot
figure; hold on, box on
imagesc(X(1,:),Y(:,1),Hfield), shading interp, view(2), colormap(gray(10))
plot(y0(:,1),y0(:,2),'-','Color',0.9*ones(1,3),'LineWidth',1);
hold on , grid on 
plot(DataStore.y1(1:end,1),DataStore.y1(1:end,2),'-','Color',[1,0,0],'LineWidth',1.2); 
plot(DataStore.y2(1:end,1),DataStore.y2(1:end,2),'-','Color',[1,1,0],'LineWidth',1.2); 
plot(DataStore.y3(1:end,1),DataStore.y3(1:end,2),'-','Color',[0,0,1],'LineWidth',1.2);
plot(DataStore.y4(1:end,1),DataStore.y4(1:end,2),'-','Color',[0,1,0],'LineWidth',1.2);
xlabel('x_1'), ylabel('x_2')
hold on

% Separatrices
plot(2*atan(sinh([-pi:0.1:pi])), 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) % center
plot(-2*atan(sinh([-pi:0.1:pi])), -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
plot(2*atan(sinh([-pi:0.1:pi]))-2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %left
plot(-2*atan(sinh([-pi:0.1:pi]))-2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
plot(2*atan(sinh([-pi:0.1:pi]))+2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %right
plot(-2*atan(sinh([-pi:0.1:pi]))+2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
axis([-2*pi-0.02 2*pi+0.02 -4-0.03 4+0.03])

set(gca,'xtick',[-2*pi,-pi,0,pi,2*pi],'xticklabel',{'-2\pi', '-\pi', '0', '\pi', '-2\pi'})
set(gca,'FontSize',14, 'LineWidth',1)
set(gcf,'Position',[100 100 350 250])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName,'PhasePlot','.eps']);
print('-dpng', '-loose', [path2figs,ModelName,'PhasePlot','.png']);


