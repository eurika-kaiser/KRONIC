path2data = '../Data/'; 
load([path2data,[ModelName1,'Ensemble.mat']])
path2figs = '../Figures/HAMILTONIAN_DUFFING/'; mkdir(path2figs)

%% Show Phase plot
[X,Y] = meshgrid([-2.4:0.01:2.4], [-2.4:0.01:2.4]);
Hfield = (1/2)*Y.^2-(X.^2)/2 + (V/4)*X.^4;% 

cmap = gray(2*11);
cmap = cmap(1:2:end,:);
figure;
box on
hold on
% contourf([-2.5:0.01:2.5], [-2.5:0.01:2.5],Hfield,[-1,-0.2,-0.1,0,0.1,0.5,1.5,2.5],'LineColor','none'), colormap(gray)
% sh = surf([-2.4:0.01:2.4], [-2.4:0.01:2.4],-ones(size(Hfield)),ones(size(Hfield)));shading interp,view(2)
% sh.FaceColor = [0 0 0];
sh = surf([-2.4:0.01:2.4], [-2.4:0.01:2.4],zeros(size(Hfield)));shading interp, view(2)
sh.FaceColor = [0 0 0];
contourf([-2.4:0.01:2.4], [-2.4:0.01:2.4],log(Hfield+0.2500001),[log([-0.2,-0.1,0,0.25:0.5:4]+0.25)],'LineColor','none'), colormap(cmap)
hold on

% Nx = size(yout,3);
% Ny = size(yout,2);
% for iy = 1:Ny
%     for ix = 1:Nx 
%         plot(squeeze(yout(1,iy,ix,:)),squeeze(yout(2,iy,ix,:)),'-','Color',0.7*ones(1,3),'LineWidth',2,'MarkerSize',10)
%     end
% end

Nic = size(yout_ctrl,2);
cmap = jet(Nic);
for iIC = 1:Nic
    hold on
    if any([1,2,3,6,7,8] == iIC);
%     plot3(squeeze(yout_ctrl(1,iIC,:)),squeeze(yout_ctrl(2,iIC,:)),squeeze(ones(size(yout_ctrl(2,iIC,:)))),'-','Color','b','LineWidth',1.2)
        plot(squeeze(yout_ctrl(1,iIC,:)),squeeze(yout_ctrl(2,iIC,:)),'-','Color',cmap(3,:),'LineWidth',1.2)
    end
end

Nic = size(yout_ctrl,2);
% cmap = jet(Nic);
for iIC = 1:Nic
    hold on
    if any([1,2,4,5,7] == iIC);
        plot(squeeze(yout_ctrl_in(1,iIC,:)),squeeze(yout_ctrl_in(2,iIC,:)),'-','Color',cmap(3,:),'LineWidth',1.2)
    end
end

x1 = -sqrt(2):0.01:sqrt(2);
plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--y')
plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--y')
drawnow

xlabel('x_1'), ylabel('x_2')
set(gca,'FontSize',14, 'LineWidth',1)
axis equal
axis([-2.5 2.5 -2.5 2.5])
set(gcf,'Position',[100 100 250 250])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName1,'PhasePlot','.eps']);