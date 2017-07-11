dt = 1;

colors = [0,0,0;
    1,0,0;
    0,0,0.7;
    0.9,0.5,0.2;
    0,0.7,0];

%% Show results
LineWidth = 3;

fhandle = figure;
plot(xLQR(:,1),xLQR(:,2),'-','Color',colors(1,:),'LineWidth',LineWidth);
hold on , grid on
if all(xKOOC(1:dt:end,1)==0)==0
    plot(xKOOC(1:dt:end,1),xKOOC(1:dt:end,2),'-','Color',colors(2,:),'LineWidth',LineWidth);
end
plot(xKOOC2(1:dt:end,1),xKOOC2(1:dt:end,2),'--','Color',colors(3,:),'LineWidth',LineWidth);
if all(xFL(1:dt:end,1)==0)==0
    plot(xFL(1:dt:end,1),xFL(1:dt:end,2),':','Color',colors(4,:),'LineWidth',LineWidth);
end
plot(xNLC(1:dt:end,1),xNLC(1:dt:end,2),'-.','Color',colors(5,:),'LineWidth',LineWidth);
xlabel('x1'), ylabel('x2')
set(gca,'FontSize',16)
axis tight
ylim(ylim_vals)
xlim([-5 0.001])
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Results_PhasePlot','.eps']);
close(fhandle);

fhandle = figure;
plot(tspan,xLQR,'-','Color',colors(1,:),'LineWidth',LineWidth); hold on , grid on
if all(xKOOC(1:dt:end,1)==0)==0
    plot(tspan,xKOOC,'-','Color',colors(2,:),'LineWidth',LineWidth);
end
plot(tspan,xKOOC2,'--','Color',colors(3,:),'LineWidth',LineWidth);
if all(xFL(1:dt:end,1)==0)==0
    plot(tspan,xFL(1:dt:end,:),':','Color',colors(4,:),'LineWidth',LineWidth);
end
plot(tNLC,xNLC,'-.','Color',colors(5,:),'LineWidth',LineWidth);
xlabel('t'), ylabel('xk')
xlim([0 50])
set(gca,'xtick',[0,20,40])
set(gca,'FontSize',16)
ylim(ylim_vals)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Results_TimeSeries','.eps']);
close(fhandle);


clear ph
fhandle = figure;
ph(1) = plot(tspan,JLQR_x,'-','Color',colors(1,:),'LineWidth',LineWidth);
% ph(1) = plot(tspan,JLQR_x/10^3,'-','Color',colors(1,:),'LineWidth',LineWidth);
hold on , grid on
if all(xKOOC(1:dt:end,1)==0)==0
    ph(2) = plot(tspan,JKOOC_x,'-','Color',colors(2,:),'LineWidth',LineWidth);
end
ph(3) = plot(tspan,JKOOC2_x,'--','Color',colors(3,:),'LineWidth',LineWidth);
% ph(3) = plot(tspan,JKOOC2_x/10^3,'--','Color',colors(3,:),'LineWidth',LineWidth);
if all(xFL(1:dt:end,1)==0)==0
    ph(4) = plot(tspan,JFL_x,':','Color',colors(4,:),'LineWidth',LineWidth);
end
ph(5) = plot(tNLC,JNLC_x,'-.','Color',colors(5,:),'LineWidth',LineWidth);
% ph(5) = plot(tNLC,JNLC_x/10^3,'-.','Color',colors(5,:),'LineWidth',LineWidth);
xlabel('t'), ylabel('Jx')
axis(axis_lim)
% axis([axis_lim(1),axis_lim(2),0,15])
set(gca,'xtick',[0,20,40])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Results_Cost','.eps']);


% Legend
clear ph
fhandle = figure;
ph(1) = plot([10,20],[3,3],'-','Color',colors(1,:),'LineWidth',2);
hold on , grid on
if all(xKOOC(1:dt:end,1)==0)==0
    ph(2) = plot([10,20],[3,3],'-','Color',colors(2,:),'LineWidth',2);
end
ph(3) = plot([10,20],[3,3],'--','Color',colors(3,:),'LineWidth',2);
if all(xFL(1:dt:end,1)==0)==0
    ph(4) = plot([10,20],[3,3],':','Color',colors(4,:),'LineWidth',2);
end
ph(5) = plot([10,20],[3,3],'-.','Color',colors(5,:),'LineWidth',2);
if all(xKOOC(1:50:end,1)==0)==0
    lhandle = legend('Linearized system','KOOC', 'KRONIC', 'Feedback linearization', 'Numer. TPBV','location','none');
else
    lhandle = legend('Linearized system', 'KRONIC', 'Numer. TPBV','location','none');
end
posvals = get(lhandle,'pos');
set(lhandle,'pos',[0.418 0.5 0.2 0.1])
axis([axis_lim(1:2),0,6])
axis off
% ph(1).Color = 0.94.*ones(1,3);
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Results_Legend','.eps']);
% close(fhandle);

