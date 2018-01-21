fontsize = 15;
LineWidth = 2;

Au200b8Ai = dlmread('../data/oriAu200b8Ai.dat');
Au200b8QGP = dlmread('../data/oriAu200b8QGP.dat');
Pb2760b8Ai = dlmread('../data/oriPb2760b8Ai.dat');
Pb2760b8QGP = dlmread('../data/oriPb2760b8QGP.dat');


figure
semilogy(Au200b8QGP(:,1),Au200b8QGP(:,2),'k-',...
    Au200b8Ai(:,1),Au200b8Ai(:,2),'k--',...
    Pb2760b8QGP(:,1),Pb2760b8QGP(:,2),'k-',...
    Pb2760b8Ai(:,1),Pb2760b8Ai(:,2),'k--',...
    'LineWidth',LineWidth)

set(gca,'YMinorTick','on')
xlim([0 6.7])
ylim([10^-2 10^6])
xlabel('t (fm)','FontSize',fontsize)
%ylabel('$\mathrm{eB}\, (\mathrm{MeV}^2)$','Interpreter','latex','FontSize',fontsize)
ylabel('eB (MeV^2)','FontSize',fontsize)
legend({'QGP Response','Vacuum'}, 'fontsize',fontsize+3)
text(5.05,100,'$\sqrt{S}=200\,\mathrm{GeV}$','Interpreter','latex','fontsize',fontsize-1)
text(5.05,5.76,'$\sqrt{S}=2760\,\mathrm{GeV}$','Interpreter','latex','fontsize',fontsize-1)
text(5.05,1.14,'$\sqrt{S}=200\,\mathrm{GeV}$','Interpreter','latex','fontsize',fontsize-1)
text(5.05,0.2727,'$\sqrt{S}=2760\,\mathrm{GeV}$','Interpreter','latex','fontsize',fontsize-1)
set(gca,'linewidth',LineWidth) 
set(gca,'fontsize',fontsize+2)