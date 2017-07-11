figure (1)
ax1 = subplot('position', [0.08 0.41 0.15 0.5]);
plot(y0(:,1),y0(:,2),'-','Color',0.7*ones(1,3),'LineWidth',1.2);
hold on , grid on 
plot(y1(:,1),y1(1:end,2),'b-','LineWidth',1.2); 
xlabel('x_1'), ylabel('x_2')
set(gca,'FontSize',14)
axis([-4 4 -2.5 2.5])
% set(gca,'xtick',[-0.4:0.4:0.4])


valmin = min([min(y0,[],1),min(y1,[],1)]);
valmax = max([max(y0,[],1),max(y1,[],1)]);
% tmp = [mod(y1(:,1),2*pi),y1(:,2)];
% valmin = min([min(y0,[],1),min(tmp,[],1)]);
% valmax = max([max(y0,[],1),max(tmp,[],1)]);
dval = 0.1*(valmax-valmin);
ax2 = subplot('position', [0.31 0.41 0.15 0.5]); 
plot(tspan,y0,'-','Color',0.7*ones(1,3),'LineWidth',1.2); hold on , grid on 
plot(tspan,y1(:,1),'b-','LineWidth',1.2); 
% plot(tspan,tmp(:,1),'b-','LineWidth',1.2); 
plot(tspan,y1(:,2),'b-','LineWidth',1.2); 
xlabel('t'), ylabel('x_k')
set(gca,'FontSize',14)
axis([0 maxt valmin-dval valmax+dval])
% axis([0 maxt -4 8])
% axis([0 maxt -4 4])

ax3 = subplot('position', [0.55 0.41 0.15 0.5]);
valmin = min([min(Hvals0),min(Hvals1)]);
valmax = max([max(Hvals0),max(Hvals1)]);
dval = 0.1*(valmax-valmin);
plot(tspan,Hvals0,'-','Color',0.7*ones(1,3),'LineWidth',1.2);
hold on , grid on
plot(tspan,Hvals1,'b-','LineWidth',1.2);
xlabel('t'), ylabel('H')
lhandle = legend('Unforced dynamics', ['KOOC of Hamiltonian with H_{ref}=',num2str(REF)],'location','none');
posvals = get(lhandle,'pos');
set(lhandle,'pos',[0.4 0.06 0.2 0.1])
set(gca,'FontSize',14)
axis([0 maxt valmin-dval valmax+dval])
% axis ([0 maxt -1 1.2])%1.2
% axis ([0 maxt -1.2 -.5])

ax4 = subplot('position', [0.82 0.41 0.15 0.5]);
valmin = 0;
valmax = max(max(cumsum(Jvals1)));
dval = 0.1*(valmax-valmin);
plot(tspan,cumsum(Jvals0),'-','Color',0.7*ones(1,3),'LineWidth',1.2);
hold on , grid on
plot(tspan,cumsum(Jvals1),'-','Color',[0 0 1],'LineWidth',1.2);
xlabel('t'), ylabel('J')
set(gca,'FontSize',14)
% axis([0 maxt 0 2/3*10^4])%2, 2/3*10^4
axis([0 maxt 0 valmax+2*dval])

%linkaxes([ax2,ax3],'x')
set(gcf,'Position',[100 100 600 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_case,'Results','.eps']);

figure;
hold on, box on
plot(tspan,y1(:,1),'b-','LineWidth',1.2); %mod(y1(:,1),(2*pi)), y1(:,1)
% plot(tspan,mod(y1(:,1),(2*pi)),'b-','LineWidth',1.2);
plot(tspan,y1(:,2),'-','Color',[0 0.7 0],'LineWidth',1.2);
xlabel('t'), ylabel('x_k')
%xlim ([0 50])
legend('x_1','x_2')
set(gca,'FontSize',14)
% axis([0 200 -0.1 0.1])
set(gcf,'Position',[100 100 600 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_case,'Results_TS','.eps']);