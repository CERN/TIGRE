clear;
close all;

%% Define Geometry
geo=defaultGeometry('nVoxel',[64,64,64]','nDetector',[128,128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi-2*pi/180,180);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
projections=addCTnoise(projections);
niter=60;
%%
[imgsirt,ressirt,errorsirt]=SIRT(projections,geo,angles,niter,'groundTruth',head);
[imgcgls,rescgls,errcgls]=CGLS(projections,geo,angles,niter,'groundTruth',head);
[imglsqr,reslsqr,errlsqr]=LSQR(projections,geo,angles,niter,'groundTruth',head);
[imglsmr,reslsmr,errlsmr]=LSMR(projections,geo,angles,niter,'groundTruth',head,'lambda',0);
[imglsmr2,reslsmr2,errlsmr2]=LSMR(projections,geo,angles,niter,'groundTruth',head,'lambda',30);
[imghlsqr,reshlsqr,errhlsqr]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head);
[imgcglstv,rescglstv,errcglstv]=IRN_TV_CGLS(projections,geo,angles,niter,'groundTruth',head,'niter_outer',4);
[imgflsqr,resflsqr,errflsqr]=hybrid_fLSQR_TV(projections,geo,angles,niter,'groundTruth',head,'lambda',1);
[imgabgmres,resabgmres,errabgmres]=AB_GMRES(projections,geo,angles,niter,'groundTruth',head);
[imgbagmres,resbagmres,errbagmres]=BA_GMRES(projections,geo,angles,niter,'groundTruth',head);
[imgabgmresfdk,resabgmresfdk,errabgmresfdk]=AB_GMRES(projections,geo,angles,niter,'groundTruth',head,'backprojector','FDK');
[imgbagmresfdk,resbagmresfdk,errbagmresfdk]=BA_GMRES(projections,geo,angles,niter,'groundTruth',head,'backprojector','FDK');

%% clean up data:
errorsirt(errorsirt==0)=nan;
errcgls(errcgls==0)=nan;
errlsqr(errlsqr==0)=nan;
errlsmr(errlsmr==0)=nan;
errlsmr2(errlsmr2==0)=nan;
errhlsqr(errhlsqr==0)=nan;
errcglstv(errcglstv==0)=nan;
errflsqr(errflsqr==0)=nan;
errabgmres(errabgmres==0)=nan;
errbagmres(errbagmres==0)=nan;
errabgmresfdk(errabgmresfdk==0)=nan;
errbagmresfdk(errbagmresfdk==0)=nan;

%%
extracolors=linspecer(10);
colors=lines(10);
colors(8:10,:)=extracolors(7:9,:);
colors=(linspecer(10));
%% plots N1:
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
figure(1)
colororder(colors)

semilogy(errorsirt/norm(head(:)),'linewidth',2);
hold on; 
semilogy(errcgls/norm(head(:)),'linewidth',2);
semilogy(errlsqr/norm(head(:)),'linewidth',2);
semilogy(errlsmr/norm(head(:)),'linewidth',2);
semilogy(errlsmr2/norm(head(:)),'linewidth',2);
semilogy(errhlsqr/norm(head(:)),'linewidth',2);
% semilogy(errcglstv/norm(head(:)),'linewidth',2);
% semilogy(errflsqr/norm(head(:)),'linewidth',2);
semilogy(errabgmres/norm(head(:)),'linewidth',2);
semilogy(errbagmres/norm(head(:)),'linewidth',2);
semilogy(errabgmresfdk/norm(head(:)),'linewidth',2);
semilogy(errbagmresfdk/norm(head(:)),'linewidth',2);
xlim([1,niter]);
ylim([7*10^-2,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
title("Error norm")
legend(["SIRT","CGLS","LSQR","LSMR $\lambda=0$","LSMR $\lambda=30$","hybrid LSQR","AB-GMRES","BA-GMRES","AB-GRMES ($B_{FDK}$)","BA-GRMES ($B_{FDK}$)"],'NumColumns',2,'interpreter','latex')

% legend(["SIRT","CGLS","LSQR","LSMR $\lambda=0$","LSMR $\lambda=30$","hybrid LSQR","CGLS-TV","hybrid fLSQR TV"],'NumColumns',2,'interpreter','latex')
set(gca,'fontsize',16)
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%

figure(2)
colororder(colors)
semilogy(ressirt/norm(projections(:)),'linewidth',2);
hold on; 
semilogy(rescgls/norm(projections(:)),'linewidth',2);
semilogy(reslsqr/norm(projections(:)),'linewidth',2);
semilogy(reslsmr/norm(projections(:)),'linewidth',2);
semilogy(reslsmr2/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr/norm(projections(:)),'linewidth',2);
% semilogy(rescglstv/norm(projections(:)),'linewidth',2);
% semilogy(resflsqr/norm(projections(:)),'linewidth',2);
semilogy(resabgmres/norm(projections(:)),'linewidth',2);
semilogy(resbagmres/norm(projections(:)),'linewidth',2);
semilogy(resabgmresfdk/norm(projections(:)),'linewidth',2);
semilogy(resbagmresfdk/norm(projections(:)),'linewidth',2);
xlim([1,niter]);
ylim([9*10^-3,10^0])
xlabel("Iteration number")
ylabel("$\|Ax-b\|/\|b\|$",'interpreter','latex')
title("Residual norm")
legend(["SIRT","CGLS","LSQR","LSMR $\lambda=0$","LSMR $\lambda=30$","hybrid LSQR","AB-GMRES","BA-GMRES","AB-GRMES ($B_{FDK}$)","BA-GRMES ($B_{FDK}$)"],'NumColumns',2,'interpreter','latex')
grid on
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
set(gca,'fontsize',16)

axes('Position',[.68 .42 .2 .2])
box on

colororder(colors)
semilogy(ressirt/norm(projections(:)),'linewidth',2);
hold on; 
semilogy(rescgls/norm(projections(:)),'linewidth',2);
semilogy(reslsqr/norm(projections(:)),'linewidth',2);
semilogy(reslsmr/norm(projections(:)),'linewidth',2);
semilogy(reslsmr2/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr/norm(projections(:)),'linewidth',2);
% semilogy(rescglstv/norm(projections(:)),'linewidth',2);
% semilogy(resflsqr/norm(projections(:)),'linewidth',2);
semilogy(resabgmres/norm(projections(:)),'linewidth',2);
semilogy(resbagmres/norm(projections(:)),'linewidth',2);
semilogy(resabgmresfdk/norm(projections(:)),'linewidth',2);
semilogy(resbagmresfdk/norm(projections(:)),'linewidth',2);
xlim([22,32]);
ylim([1.98*10^-2,2.05*10^-2])
set(gca,'Xtick',[22,26,30,32])
grid on
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
set(gca,'fontsize',10)

%%
set(0,'defaultTextInterpreter','none');

plotImg(cat(3,[imghlsqr,imgabgmres,imgbagmres,imgabgmresfdk,imgbagmresfdk],[imgsirt,imgcgls,imglsqr,imglsmr,imglsmr2]),'dim',1,'clims',[0 1],'slice',32)
plotImg(cat(3,[head-imghlsqr,head-imgabgmres,head-imgbagmres,head-imgabgmresfdk,head-imgbagmresfdk],[head-imgsirt,head-imgcgls,head-imglsqr,head-imglsmr,head-imglsmr2]),'dim',1,'clims',[-0.1 0.1],'slice',32,'colormap',redblue)

% plotImg([head-imgsirt,head-imgcgls,head-imglsqr,head-imglsmr,head-imglsmr2,head-imghlsqr,head-imgcglstv,head-imgflsqr],'dim',1,'clims',[-0.1 0.1],'slice',32,'colormap',redblue)

%%
geo=defaultGeometry('nVoxel',[64,64,64]','nDetector',[128,128]);                     
nangles=60;
angles=linspace(0,2*pi-2*pi/nangles,nangles);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
projections=addCTnoise(projections,'Poisson',1e5);
niter=60;
%% Experiment 2
[imglsqr,reslsqr,errlsqr]=LSQR(projections,geo,angles,niter,'groundTruth',head);
%GCV
[imghlsqr1,reshlsqr1,errhlsqr1]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head);
%Discrepacy
[imghlsqr2,reshlsqr2,errhlsqr2]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head,'NoiseLevel',0.015);

% lambda
[imghlsqr3,reshlsqr3,errhlsqr3]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head,'lambda',20);
[imghlsqr4,reshlsqr4,errhlsqr4]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head,'lambda',2);
[imghlsqr5,reshlsqr5,errhlsqr5]=hybrid_LSQR(projections,geo,angles,niter,'groundTruth',head,'lambda',200);

%%
errlsqr(errlsqr==0)=nan;
errhlsqr1(errhlsqr1==0)=nan;
errhlsqr2(errhlsqr2==0)=nan;
errhlsqr3(errhlsqr3==0)=nan;
errhlsqr4(errhlsqr4==0)=nan;
errhlsqr5(errhlsqr5==0)=nan;

%% plots N1:
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
figure(1)
colororder(lines(6))

semilogy(errlsqr/norm(head(:)),'linewidth',2);
hold on; 
semilogy(errhlsqr1/norm(head(:)),'linewidth',2);
semilogy(errhlsqr2/norm(head(:)),'linewidth',2);
semilogy(errhlsqr3/norm(head(:)),'linewidth',2);
semilogy(errhlsqr4/norm(head(:)),'linewidth',2);
semilogy(errhlsqr5/norm(head(:)),'linewidth',2);

xlim([1,niter]);
ylim([7*10^-2,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
title("Error norm")
legend(["LSQR","hybrid LSQR (GCV)","hybrid LSQR (DP)","hybrid LSQR $\lambda=20$","hybrid LSQR $\lambda=2$","hybrid LSQR $\lambda=200$"],'NumColumns',2,'interpreter','latex')
set(gca,'fontsize',20)
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%
%% plots N1:
figure(2)
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
colororder(lines(6))

semilogy(reslsqr(1,:)'/norm(projections(:)),'linewidth',2);
hold on; 
semilogy(reshlsqr1(1,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr2(1,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr3(1,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr4(1,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr5(1,:)/norm(projections(:)),'linewidth',2);

xlim([1,niter]);
ylim([9*10^-3,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
ylabel("$\|Ax-b\|/\|b\|$",'interpreter','latex')
legend(["LSQR","hybrid LSQR (GCV)","hybrid LSQR (DP)","hybrid LSQR $\lambda=20$","hybrid LSQR $\lambda=2$","hybrid LSQR $\lambda=200$"],'NumColumns',2,'interpreter','latex')
set(gca,'fontsize',20)
title("Residual norm (implicit)")

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%
figure(3)
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
colororder(lines(5))

semilogy(reslsqr(2,:)'/norm(projections(:)),'linewidth',2);
hold on; 
semilogy(reshlsqr1(2,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr2(2,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr3(2,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr4(2,:)/norm(projections(:)),'linewidth',2);
semilogy(reshlsqr5(2,:)/norm(projections(:)),'linewidth',2);

xlim([1,niter]);
ylim([9*10^-3,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
ylabel("$\|Ax-b\|/\|b\|$",'interpreter','latex')
legend(["LSQR","hybrid LSQR (GCV)","hybrid LSQR (DP)","hybrid LSQR $\lambda=20$","hybrid LSQR $\lambda=2$","hybrid LSQR $\lambda=200$"],'NumColumns',2,'interpreter','latex')
set(gca,'fontsize',20)
title("Residual norm (explicit)")

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%
plotImg([imglsqr,imghlsqr1,imghlsqr2,imghlsqr3,imghlsqr4,imghlsqr5],'dim',1,'clims',[0 1],'slice',32)
plotImg([head-imglsqr,head-imghlsqr1,head-imghlsqr2,head-imghlsqr3,head-imghlsqr4,head-imghlsqr5],'dim',1,'clims',[-0.1 0.1],'slice',32,'colormap',redblue)


%% TV EXPERIMENTS
geo=defaultGeometry('nVoxel',[64,64,64]','nDetector',[128,128]);                     
nangles=60;
angles=linspace(0,2*pi-2*pi/nangles,nangles);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
projections=addCTnoise(projections,'Poisson',1e4);

%%
niter=60;
[imgcgls,rescgls,errcgls]=CGLS(projections,geo,angles,niter,'groundTruth',head);
[imgcglstv,rescglstv,errcglstv]=IRN_TV_CGLS(projections,geo,angles,niter,'groundTruth',head,'niter_outer',4);
[imgflsqr,resflsqr,errflsqr]=hybrid_fLSQR_TV(projections,geo,angles,niter,'groundTruth',head,'lambda',5);
%%
errcgls(errcgls==0)=nan;
errcglstv(errcglstv==0)=nan;
errflsqr(errflsqr==0)=nan;


%% plots N1:
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
figure(1)
colororder(lines(6))

semilogy(errcgls/norm(head(:)),'linewidth',2);
hold on; 
semilogy(errcglstv/norm(head(:)),'linewidth',2);
semilogy(errflsqr/norm(head(:)),'linewidth',2);


xlim([1,niter]);
ylim([7*10^-2,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
title("Error norm")
legend(["CGLS","CGLS-TV","hybrid fLSQR TV"],'interpreter','latex')
set(gca,'fontsize',16)
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%
%% plots N1:
figure(2)
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Helvetica')
colororder(lines(6))

semilogy(rescgls/norm(projections(:)),'linewidth',2);
hold on; 
semilogy(rescglstv/norm(projections(:)),'linewidth',2);
semilogy(resflsqr/norm(projections(:)),'linewidth',2);

xlim([1,niter]);
ylim([9*10^-3,10^0])
xlabel("Iteration number")
ylabel("$\|x-x_{gt}\|/\|x_{gt}\|$",'interpreter','latex')
ylabel("$\|Ax-b\|/\|b\|$",'interpreter','latex')
legend(["CGLS","CGLS-TV","hybrid fLSQR TV"],'interpreter','latex')
set(gca,'fontsize',16)
title("Residual norm")

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9/1.2, 7/1.2], 'PaperUnits', 'Inches', 'PaperSize', [9/1.2 7/1.2])
grid on
%%
plotImg(cat(3,[imgcgls,imgcglstv,imgflsqr]),'dim',1,'clims',[0 1],'slice',32)
plotImg([head-imgcgls,head-imgcglstv,head-imgflsqr],'dim',1,'clims',[-0.1 0.1],'slice',32,'colormap',redblue)