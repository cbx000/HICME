%% read data
Au20001 = csvread('../data/ratioAu200GeV0.1.dat',1,0);
Au20002 = csvread('../data/ratioAu200GeV0.2.dat',1,0);
Au20003 = csvread('../data/ratioAu200GeV0.3.dat',1,0);

%% plot
linewidth = 2;
fontsize = 18;
markersize = 10;

figure
box on
plot(Au20001(:,1),Au20001(:,2),'-.',...
     Au20002(:,1),Au20002(:,2),'--',...
     Au20003(:,1),Au20003(:,2),'-','LineWidth',linewidth);
set(gca,'linewidth',2,'FontName','Times','FontSize',fontsize-2);
legend({'$\lambda = 0.1 R$','$\lambda = 0.2 R$', '$\lambda = 0.3 R$'},...
     'Interpreter','latex','Location','northwest');
xlabel('$b/R$','Interpreter','latex')
ylabel('$|a_{+-}|/a_{++}$','Interpreter','latex')
set(gca,'FontName','Times','FontSize',fontsize-2)


