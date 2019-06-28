load sentry498_sensors
for ii = 0:17421
    jj = ii + 6748;
    kk = ii + 12447;
    ll = ii+1;
    ssp(ll,1) = depth.depth(1,jj);
    ssp(ll,2) = svp.sound_speed(1,kk);
end
nn = 0;
oo = 1;
for mm = 1:length(ssp(:,1))
    if ssp(mm,1) > nn
        cmat(oo,1) = round(ssp(mm,1));
        cmat(oo,2) = ssp(mm,2);
        nn = nn+1;
        oo = oo+1;
    else
    end
end

% Sound Speed Figure
soundspeedFig=figure;
figure(soundspeedFig)
plot(cmat(:,2),cmat(:,1))
set(gca,'ydir','reverse')
xlabel('m/s'); ylabel('depth (m)');