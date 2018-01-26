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

Pb276001 = csvread('../data/Pb2760GeV0.1.dat',1,0);
Pb276002 = csvread('../data/Pb2760GeV0.2.dat',1,0);
Pb276003 = csvread('../data/Pb2760GeV0.3.dat',1,0);

Au20001Diff = Au20001(2:8,1) - Au20001(2:8,2);
Au20002Diff = Au20002(2:8,1) - Au20002(2:8,2);
Au20003Diff = Au20003(2:8,1) - Au20003(2:8,2);

Pb276001Diff = Pb276001(1:8,1) - Pb276001(1:8,2);
Pb276002Diff = Pb276002(1:8,1) - Pb276002(1:8,2);
Pb276003Diff = Pb276003(1:8,1) - Pb276003(1:8,2);

linewidth = 2;
fontsize = 18;
markersize = 10;
%% plot exp H

figure
hold on
box on
markdersize = 9;
plot(2:8,Au200Hsame,'bo','MarkerFaceColor','b','MarkerSize',markdersize)
plot(2:8,Au200Hopp,'bo','MarkerSize',markdersize)
plot(1:8,Pb2760Hsame,'rs','MarkerFaceColor','r','MarkerSize',markdersize)
plot(1:8,Pb2760Hopp,'rs','MarkerSize',markdersize)
set(gca,'linewidth',2);
legend({'Au $200\,\mathrm{GeV}\,H_{\mathrm{SS}}$','Au $200\,\mathrm{GeV}\,H_{\mathrm{OS}}$',...
    'Pb $2760\,\mathrm{GeV}\,H_{\mathrm{SS}}$','Pb $2760\,\mathrm{GeV}\,H_{\mathrm{OS}}$'},...
    'Interpreter','latex','Location','northwest')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',fontsize-2)
xlim([0.5 8.5])
set(gca,'XTickLabel',{'0-5%','5-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%'})
xlabel('Centrality','FontSize',fontsize)
ylabel('$H_{\mathrm{SS}}$ or $H_{\mathrm{OS}}$','Interpreter','latex','FontSize',fontsize)

%% compare to Hdiff at RHIC
Au200HAlpha = 379.598509214; % lambda = 0.3
Pb2760HAlpha = 55490.6144561 ; % lambda = 0.3
% linewidth = 2;
% fontsize = 18;
% markersize = 10;
figure
hold on
box on
plot(2:8,Au200Hdiff,'ro','MarkerFaceColor','r','MarkerSize',markersize)
plot(1:8,(Au20003(1:8,1) - Au20003(1:8,2))*Au200HAlpha,'-k','LineWidth',linewidth)
set(gca,'linewidth',2);
legend({'Au $200\,\mathrm{GeV}$ Exp.','Au $200\,\mathrm{GeV}$ Theory'},...
    'Interpreter','latex','Location','northwest')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',fontsize-2)
xlim([0.5 8.5])
set(gca,'XTickLabel',{'0-5%','5-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%'})
xlabel('Centrality','FontSize',fontsize)
ylabel('$H_{\mathrm{SS}}-H_{\mathrm{OS}}$','Interpreter','latex','FontSize',fontsize)




%% compare to Hdiff at LHC
linewidth = 2;
fontsize = 18;
markersize = 10;
figure
hold on
box on
plot(1:8,Pb2760Hdiff,'ro','MarkerFaceColor','r','MarkerSize',markersize)
plot(1:8,Pb276003Diff*Pb2760HAlpha,'-k','LineWidth',linewidth)
set(gca,'linewidth',2);
legend({'Pb $2760\,\mathrm{GeV}$ Exp.','Pb $2760\,\mathrm{GeV}$ Theory'},...
    'Interpreter','latex','Location','northwest')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',fontsize-2)
xlim([0.5 8.5])
set(gca,'XTickLabel',{'0-5%','5-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%'})
xlabel('Centrality','FontSize',fontsize)
ylabel('$H_{\mathrm{SS}}-H_{\mathrm{OS}}$','Interpreter','latex','FontSize',fontsize)

%% plot all
linewidth = 2;
fontsize = 18;
markersize = 10;
figure
hold on
box on
plot(2:7,Au200Hdiff(1:end-1),'ro','MarkerFaceColor','r','MarkerSize',markersize)
plot(1:7,(Au20003(1:7,1) - Au20003(1:7,2))*Au200HAlpha,'-r','LineWidth',linewidth)
plot(1:7,Pb2760Hdiff(1:end-1),'bo','MarkerSize',markersize)
plot(1:7,Pb276003Diff(1:end-1)*Pb2760HAlpha,'--b','LineWidth',linewidth)
set(gca,'linewidth',2);
legend({'Au-Au $200\,\mathrm{GeV}$ Exp.','Au-Au $200\,\mathrm{GeV}$ Theory',...
    'Pb-Pb $2760\,\mathrm{GeV}$ Exp.','Pb-Pb $2760\,\mathrm{GeV}$ Theory'},...
    'Interpreter','latex','Location','northwest')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',fontsize-2)
xlim([0.5 7.5])
set(gca,'XTickLabel',{'0-5%','5-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%'})
xlabel('Centrality','FontSize',fontsize)
ylabel('$H_{\mathrm{SS}}-H_{\mathrm{OS}}$','Interpreter','latex','FontSize',fontsize)

