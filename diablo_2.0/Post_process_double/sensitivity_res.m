
a = load('PSD_tke_Pr9_Ri16_D_3_res1.mat');
b = load('PSD_tke_Pr9_Ri16_D_3_res3.mat');

 figure('position',[1000,1000,750,600]); 
t= tiledlayout(2,2,"TileSpacing",'compact','Padding','compact');

nexttile;
loglog(a.ks,a.TKE(:,1).*a.ks'.^3,'linewidth',2); hold on;
loglog(b.ks,b.TKE(:,1).*b.ks'.^3,'linewidth',2);
%hold on; plot(a.ks,a.ks.^(-5/3),'linewidth',1.5,'color','k')
set(gca,'fontsize',14,'linewidth',1.2,'xtick',[1e0,1e1,1e2]);
%xlabel('k');
ylabel('k^3E_{uu}(k)');
ylim([1e-10,1e1])
title('t=t_{2d}','fontweight','normal','fontangle','italic');

nexttile;
loglog(a.ks,a.TKE(:,2).*a.ks'.^3,'linewidth',2); hold on; 
loglog(b.ks,b.TKE(:,2).*b.ks'.^3,'linewidth',2);
%hold on; plot(a.ks,a.ks.^(-5/3),'linewidth',1.5,'color','k')
set(gca,'fontsize',14,'linewidth',1.2,'xtick',[1e0,1e1,1e2]);
ylim([1e-10,1e1])
title('t=t_{2d}+24','fontweight','normal','fontangle','italic');

nexttile;
loglog(a.ks,a.Sp_kb(:,1).*a.ks'.^3,'linewidth',2); hold on;
loglog(b.ks,b.Sp_kb(:,1).*b.ks'.^3,'linewidth',2);
%hold on; plot(a.ks,a.ks.^(-5/3),'linewidth',1.5,'color','k')
set(gca,'fontsize',14,'linewidth',1.2,'xtick',[1e0,1e1,1e2]);
xlabel('k');
ylabel('k^3E_{bb}(k)');
ylim([1e-5,1e1])

nexttile;
loglog(a.ks,a.Sp_kb(:,2).*a.ks'.^3,'linewidth',2); hold on;
loglog(b.ks,b.Sp_kb(:,2).*b.ks'.^3,'linewidth',2);
%hold on; plot(a.ks,a.ks.^(-5/3),'linewidth',1.5,'color','k')
set(gca,'fontsize',14,'linewidth',1.2,'xtick',[1e0,1e1,1e2]);
xlabel('k');%ylabel('E_{bb}(k)');
ylim([1e-5,1e1]);
legend('\Delta x = 0.018','\Delta x = 0.0285 (Smyth & Moum 2000)','fontsize',12,'location','southwest','box','off');
