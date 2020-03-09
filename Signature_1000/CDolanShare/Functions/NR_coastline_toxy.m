% NR coastline to xy

SL = load('c:\Users\wmkra\Dropbox\north_river\north_river_zodiac\li\shoreline_northriver.mat');

close all
figure
subplot(2,2,1);

plot(SL.slnrx,SL.slnry,'k-'); hold on

%find selection around bend
SL.part1 = [5733:5803];
SL.part2 = [3010:3032];
SL.part3 = [3161:3200];
SL.part = [SL.part1,SL.part2,SL.part3];

plot(SL.slnrx(SL.part1),SL.slnry(SL.part1),'r-');
plot(SL.slnrx(SL.part2),SL.slnry(SL.part2),'c-');
plot(SL.slnrx(SL.part3),SL.slnry(SL.part3),'b-');

xlim([-70.77 -70.75]);

% translate to rivercoordinates
[SL.xNR,SL.yNR] = northriver_xy_WK(SL.slnry,SL.slnrx); %x = -70.... = oosterlengte

subplot(2,2,2);
% plot(SL.xNR,SL.yNR,'k-'); hold on
plot(SL.xNR(SL.part1),SL.yNR(SL.part1),'r.-'); hold on
plot(SL.xNR(SL.part2),SL.yNR(SL.part2),'c.-');
plot(SL.xNR(SL.part3),SL.yNR(SL.part3),'b.-');

SL.xNRsouth = SL.xNR([SL.part1]);
SL.yNRsouth = SL.yNR([SL.part1]);
SL.xNRnorth = SL.xNR([SL.part2,SL.part3]);
SL.yNRnorth = SL.yNR([SL.part2,SL.part3]);

save('Coastline_xNR_yNR.mat','SL');