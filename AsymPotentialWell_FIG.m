path2data = '../Data/'; ModelName1 = ['AsymPotentialWell_','B11_'];
load([path2data,[ModelName1,'Data.mat']])
path2figs = '../Figures/ASYM_POT_WELL/'; mkdir(path2figs)


%% Show Potential function
V = @(x)(x.^4/4 - x.^2/2 - (a/3).*x.^3 + a*x);
figure,hold on, box on
x = [-3:0.01:3]; 
plot(x,x.^4/4 - x.^2/2 - (a/3).*x.^3 + a*x,'-k','LineWidth',3);
xlim([-2.5 2.5]),ylim([-0.5 1])
plot(x0(1,1),V(x0(1,1)),'o','Color',[0,0.5,1],'MarkerFaceColor',[0,0.5,1]);
plot(x0(2,1),V(x0(2,1)),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1]);
plot(-1,V(-1),'or','MarkerFaceColor','r');
plot(x,V(a).*ones(length(x)),'--k')
plot(a,V(a),'or','MarkerFaceColor','r');
plot(1,V(1),'or','MarkerFaceColor','r');
xlabel('x1'), ylabel('V')
set(gca,'ytick',[0,1])
set(gca,'FontSize',14, 'LineWidth',1)
set(gcf,'Position',[100 100 250 120])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'Potential','.eps']);


%% Show phase plot with controlled trajectories
cmap = gray(2*11);
cmap = cmap(1:2:end,:);
LineWidth = 2;
clear ph

figure;
box on
hold on
sh = surf(x,y,-0.4.*ones(size(Hgrid)));shading interp, view(2), hold on
sh.FaceColor = [0 0 0];
contourf(x,y,log(Hgrid+0.42),[log([-0.2,-0.1,0,0.25:0.5:4]+0.42)],'LineColor','none'), colormap(cmap)
% for i = 1:3
%     ph(1) = plot3(y0(1:2000,1,i),y0(1:2000,2,i),ones(size(y0(1:2000,2,i))),'--','Color',[0.8,0.8,0],'LineWidth',1);
% end
ph(1) = plot3(y0(1:10000,1,1),y0(1:10000,2,1),ones(size(y0(1:10000,2,1))),'--','Color',[0.8,0.8,0],'LineWidth',1);
ph(1) = plot3(y0(1:5500,1,2),y0(1:5500,2,2),ones(size(y0(1:5500,2,2))),'--','Color',[0.8,0.8,0],'LineWidth',1);
ph(1) = plot3(y0(1:40000,1,3),y0(1:40000,2,3),ones(size(y0(1:40000,2,3))),'--','Color',[0.8,0.8,0],'LineWidth',1);
ph(2) = plot3(y1(:,1,1),y1(:,2,1),ones(size(y1(:,2,1))),'-','Color',[0,0.5,1],'LineWidth',LineWidth+1);
ph(3) = plot3(y1(:,1,2),y1(:,2,2),ones(size(y1(:,2,2))),'-','Color',[0,1,1],'LineWidth',1);
plot3(x0(1,1),x0(1,2),2,'o','Color',[0.8,0.8,0],'MarkerFaceColor',[0.8,0.8,0]);
% plot3(x0(2,1),x0(2,2),2,'o','Color',[0.8,0.8,0],'MarkerFaceColor',[0.8,0.8,0],'MarkerSize',4);
% plot3(xREF(1),xREF(2),2,'ok','MarkerFaceColor','k');
plot3(-1,0,2,'ok','MarkerFaceColor','r');
plot3(a,0,2,'ok','MarkerFaceColor','r');
plot3(1,0,2,'ok','MarkerFaceColor','r');
axis equal
plot([-2.5 -2.5],[-1.5 1.5],'-k')
plot([-2.5 2.5],[1.5 1.5],'-k')
xlim([-2.5 2.5]), ylim([-1.5 1.5])
% legend(ph,'Unforced','Controlled')
xlabel('x1'), ylabel('x2')
set(gca,'FontSize',14, 'LineWidth',1)
set(gcf,'Position',[100 100 250 180])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName1,'PhasePlot','.eps']);

%% Show cost, actuation etc. 
LineWidth = 2;
fhandle = figure; hold on, box on, grid on
plot(tspan,Hvals0,'-','Color',[0.8,0.8,0],'LineWidth',2);
plot(tspan,Hvals1(:,1),'-','Color',[.1,.1,1],'LineWidth',2)
plot(tspan,Hvals1(:,2),'--','Color',[0,1,1],'LineWidth',2)
xlabel('t'), ylabel('H'), ylim([-0.41 0.81]), set(gca,'ytick',[-0.4:0.4:0.8])
set(gca,'FontSize',16)
axis tight
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'H','.eps']);


figure; hold on, box on, grid on
plot(tspan,Jvals0(:,1),'-','Color',[0.8,0.8,0],'LineWidth',2);
plot(tspan,Jvals1(:,1),'-','Color',[.1,.1,1],'LineWidth',2);
plot(tspan,Jvals1(:,2),'--','Color',[0,1,1],'LineWidth',2);
xlabel('t'), ylabel('J')
set(gca,'FontSize',16)
axis tight
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'J','.eps']);

figure; hold on
plot([0,1],[0,0],'-','Color',[0.8,0.8,0],'LineWidth',6);
plot([0,1],[0,0],'-','Color',[.1,.1,1],'LineWidth',6);
plot([0,1],[0,0],'-','Color',[0,1,1],'LineWidth',6);
pl = legend('Unforced', 'Controlled1','Controlled2'); pl.Position = [0.3 0.35 0.5 0.2];
xlim([-1 2]), axis off
set(gca,'FontSize',16)
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'legend','.eps']);

clear ph
fhandle = figure; hold on, box on, grid on
% plot(tspan,uvals0,'-','Color',[0.8,0.8,0],'LineWidth',2);
plot(tspan,uvals1(1,:,1),'-','Color',[.1,.1,1],'LineWidth',2);
plot(tspan,uvals1(2,:,1),'-','Color',[.1,0.7,.1],'LineWidth',2);
xlabel('t'), ylabel('u'),ylim([-2.1 2.1]), xlim([0 50])
set(gca,'FontSize',16)
legend('u1','u2')
set(gca,'xtick',[0,20,40])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName1,'_u','.eps']);



