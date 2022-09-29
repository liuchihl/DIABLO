cc=1:2:300;
cl=jet(length(cc));
% close all;

figure; 

for i=1:length(cc)
    plot(uv(:,i),gyf,'color',cl(i,:)); hold on;
end
% hold on; figure;plot(-sech(gyf).^2,gyf,'k','linewidth',1.4)

%%
figure; [cc,hh]=contourf(tii,gyf,uv,100);set(hh,'edgecolor','none')
xlabel('t'); ylabel('$\left<u''w''\right>_{xy}$','interpreter','latex');