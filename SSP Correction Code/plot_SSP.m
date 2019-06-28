function plot_SSP(cmat)
% Sound Speed Figure
soundspeedFig=figure;
figure(soundspeedFig)
plot(cmat(:,2),cmat(:,1))
set(gca,'ydir','reverse')
xlabel('m/s'); ylabel('depth (m)');
end