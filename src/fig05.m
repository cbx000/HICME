
% H data
Au200Hsame = csvread('../data/Au200Hsame.txt');
Au200Hopp = csvread('../data/Au200Hopp.txt');
Pb2760Hsame = csvread('../data/Pb2760Hsame.txt');
Pb2760Hopp = csvread('../data/Pb2760Hopp.txt');

Au200Hdiff = Au200Hsame - Au200Hopp;
Pb2760Hdiff = Pb2760Hsame - Pb2760Hopp;

% theory data
Au20001 = csvread('../data/Au200GeV0.1.dat',1,0);
Au20002 = csvread('../data/Au200GeV0.2.dat',1,0);
Au20003 = csvread('../data/Au200GeV0.3.dat',1,0);
Cu20001 = csvread('../data/Cu200GeV0.1.dat',1,0);
Cu20002 = csvread('../data/Cu200GeV0.2.dat',1,0);
Cu20003 = csvread('../data/Cu200GeV0.3.dat',1,0);


%% compare Au Cu Hdiff at RHIC
Au200HAlpha = 379.598509214; % lambda = 0.3
figure
hold on
box on
plot(1:7,(Au20003(1:7,1) - Au20003(1:7,2))*Au200HAlpha,'-or','LineWidth',linewidth)
plot(1:7,(Cu20003(1:7,1) - Cu20003(1:7,2))*Au200HAlpha,'-ob','LineWidth',linewidth)
set(gca,'linewidth',2);
legend({'Au-Au $200\,\mathrm{GeV}$ Theory','Cu-Cu $200\,\mathrm{GeV}$ Theory'},...
    'Interpreter','latex','Location','northwest')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',fontsize-2)
xlim([0.5 7.5])
set(gca,'XTickLabel',{'0-5%','5-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%'})
xlabel('Centrality','FontSize',fontsize)
ylabel('$a_{++}-a_{+-}$','Interpreter','latex','FontSize',fontsize)
