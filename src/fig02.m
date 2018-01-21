fontsize = 15;
LineWidth = 2;

skokov{1} = csvread('../data/skokov2014fig1One.txt');
skokov{2} = csvread('../data/skokov2014fig1Two.txt');
skokov{3} = csvread('../data/skokov2014fig1Three.txt');
skokov{4} = csvread('../data/skokov2014fig1Four.txt');

tuchin.blue = csvread('../data/tuchin2013fig4blue.txt');
tuchin.brown = csvread('../data/tuchin2013fig4brown.txt');
tuchin.green = csvread('../data/tuchin2013fig4green.txt');
tuchin.red = csvread('../data/tuchin2013fig4red.txt');
tuchin.black2016 = csvread('../data/tuchin2016fig2.txt');


Au200b8Ai = dlmread('../data/oriAu200b8Ai.dat');
Au200b8QGP = dlmread('../data/oriAu200b8QGP.dat');
Au200b7QGP = dlmread('../data/oriAu200b7QGP.dat');
Au200b6QGP = dlmread('../data/oriAu200b6QGP.dat');
Pb2760b8Ai = dlmread('../data/oriPb2760b8Ai.dat');
Pb2760b8QGP = dlmread('../data/oriPb2760b8QGP.dat');


hbar = 197.32696;
mpi = 140.0;
tuchin.blue(:,2) = tuchin.blue(:,2)*hbar^2;
tuchin.brown(:,2) = tuchin.brown(:,2)*hbar^2;
tuchin.green(:,2) = tuchin.green(:,2)*hbar^2;
tuchin.red(:,2) = tuchin.red(:,2)*hbar^2;
tuchin.black2016(:,2) = tuchin.black2016(:,2)*hbar^2;

skokov{1}(:,2) = skokov{1}(:,2)*mpi^2;
skokov{2}(:,2) = skokov{2}(:,2)*mpi^2;
skokov{3}(:,2) = skokov{3}(:,2)*mpi^2;
skokov{4}(:,2) = skokov{4}(:,2)*mpi^2;

% RHIC ???? b=6 ?????????? skokov????
figure
subplot(2,1,1)
p=semilogy(Au200b6QGP(1:251,1),Au200b6QGP(1:251,2),'k-',...
    skokov{1}(25:end,1),skokov{1}(25:end,2),'k-.',...
    skokov{2}(26:end,1),skokov{2}(26:end,2),'k:',...
    'LineWidth',LineWidth);
text(0.1,3e5,'(a)','FontSize',fontsize)
xlabel('t (fm)','FontSize',fontsize)
ylabel('eB (MeV^2)','FontSize',fontsize)
h = legend('Our method', 'Skokov: in the vacuum','Skokov: $\sigma_{\mathrm{LQCD}}$');
set(h,'Interpreter','latex')
set(gca,'linewidth',LineWidth) 
set(gca,'fontsize',fontsize)
% RHIC ???? b=7 ?????????? Tuchin????
subplot(2,1,2)
p2=semilogy(Au200b7QGP(:,1),Au200b7QGP(:,2),'k-',...
    tuchin.blue(1:53,1),tuchin.blue(1:53,2),'k--',...
    tuchin.red(1:141,1),tuchin.red(1:141,2),'k:',...
    tuchin.black2016(1:105,1),tuchin.black2016(1:105,2),'k-.',...
    'LineWidth',LineWidth);
text(0.2,3e5,'(b)','FontSize',fontsize)
xlabel('t (fm)','FontSize',fontsize)
ylabel('eB (MeV^2)','FontSize',fontsize)
h = legend('Our method', 'Tuchin: in the vacuum','Tuchin: $\sigma = 5.8\,\mathrm{MeV}$',...
    'Tuchin: $B_{\mathrm{init}} + B_{\mathrm{val}}$');
set(h,'Interpreter','latex')
set(gca,'linewidth',LineWidth) 
set(gca,'fontsize',fontsize)