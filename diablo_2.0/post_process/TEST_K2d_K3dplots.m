figure; 
semilogy(time,mean(k2d),'b','linewidth',1.5);
hold on;
semilogy(tt{1},K2d{1},'ob','linewidth',1.5);
semilogy(time,mean(k3d),'r','linewidth',1.5);
semilogy(tt{1},K3d{1},'or','linewidth',1.5);
semilogy(time,mean(tke),'--k','linewidth',1.5);

legend('$K_{2d}$ (mean.h5)','$K_{2d}$ (3D flow)','$K_{3d}$ (mean.h5)',...
    '$K_{3d}$ (3D flow)','TKE','fontsize',16,'interpreter','latex','location','southeast');
xlabel('$t$','fontsize',16,'interpreter','latex')
title('$Pr=1, Ri_0=0.1, Re=200, d=10, KICK=0.1$','interpreter','latex');
set(gca,'fontsize',12,'xminortick','on','yminortick','on');

