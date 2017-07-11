fhandle = figure;
box on
hold on
for iy = 1:1:Ny0
    for ix = 1:1:Nx0 
        plot(squeeze(yout(1,iy,ix,:)),squeeze(yout(2,iy,ix,:)),'-','Color',0.7*ones(1,3),'LineWidth',2,'MarkerSize',10)
    end
end
plot(x,(1/slope_stab_man)*x.^2,'--k','LineWidth',2)

for iy = 1:1:Ny
    for ix = 1:1:Nx 
        [t,xKOOC2] = ode45(vf3,[0:0.01:60],squeeze(yIC(:,iy,ix))');
        plot(xKOOC2(:,1),xKOOC2(:,2),'-','Color',[0 0 0.5],'LineWidth',2,'MarkerSize',10)
    end
end    

axis ([-4 4 -5 20])
xlabel('x1')
ylabel('x2')
set(gca,'FontSize',16)

%linkaxes([ax2,ax3],'x')
set(gcf,'Position',[100 100 250 200])
set(gcf,'PaperPositionMode','auto')
print('-painters', '-depsc2', '-loose', [ModelName,'PhasePlot','.eps']);

close(fhandle);