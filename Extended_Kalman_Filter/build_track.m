function [track]=build_track(freq,steps)

% freq=10;
% steps=4000*freq;

track.stw=zeros(steps,1);
track.hdg=zeros(steps,1);

%% Straight line
% track.stw(:)=.5;
% track.hdg(:)=30;
%% Multiple 90 degree turns
% for ii=1:steps
%     track(ii,1)=.5;
% %Define courses for simulation
%     if ii<=400*freq
%         track(ii,2)=30;
%     elseif ii<=490*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii<=800*freq
%         track(ii,2)=120;
%     elseif ii<=860*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii==860*freq+1
%         track(ii,2)=track(ii-1,2)*-1+.1;
%     elseif ii<=890*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii<=1200*freq
%         track(ii,2)=-150;
%     elseif ii<=1230*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     elseif ii==1230*freq+1
%         track(ii,2)=track(ii-1,2)*-1-.1;
%     elseif ii<=1290*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     elseif ii<=1600*freq
%         track(ii,2)=120;
%     elseif ii<=1690*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     elseif ii<=2400*freq
%         track(ii,2)=30;
%     elseif ii<=2490*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii<=2800*freq
%         track(ii,2)=120;
%     elseif ii<=2860*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii==2860*freq+1
%         track(ii,2)=track(ii-1,2)*-1+.1;
%     elseif ii<=2890*freq
%         track(ii,2)=track(ii-1,2)+.1;
%     elseif ii<=3200*freq
%         track(ii,2)=-150;
%     elseif ii<=3230*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     elseif ii==3230*freq+1
%         track(ii,2)=track(ii-1,2)*-1-.1;
%     elseif ii<=3290*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     elseif ii<=3600*freq
%         track(ii,2)=120;
%     elseif ii<=3690*freq
%         track(ii,2)=track(ii-1,2)-.1;
%     else
%         track(ii,2)=30;
%     end
%% Racetrack for bias determination
% for ii=1:steps
%     if ii<=steps*.33
%         track.stw(ii)=.5;
%         track.hdg(ii)=30;
%     elseif ii<=(steps*.33)+(180*freq)
%         track.stw(ii)=.5;
%         track.hdg(ii)=track.hdg(ii-1)-.1;
%     elseif ii<=steps*.66
%         track.stw(ii)=.5;
%         track.hdg(ii)=-150;
%     elseif ii<=(steps*.66)+(180*freq)
%         track.stw(ii)=.5;
%         track.hdg(ii)=track.hdg(ii-1)-.1;
%     else 
%         track.stw(ii)=.5;
%         track.hdg(ii)=30;
%     end
% end
%     track.hdg=90.-track.hdg;
%     for jj=1:length(track.hdg)
%         if track.hdg(jj)>180
%             track.hdg(jj)=track.hdg(jj)-360;
%         elseif track.hdg(jj)<-180
%             track.hdg(jj)=track.hdg(jj)+360;
%         end
%     end
%% Out and back
for ii=1:steps
    if ii<=steps*.5
        track.stw(ii)=.5;
        track.hdg(ii)=30;
    elseif ii<=(steps*.5)+(180*freq)
        track.stw(ii)=.5;
        track.hdg(ii)=track.hdg(ii-1)-(1/freq);
    else 
        track.stw(ii)=.5;
        track.hdg(ii)=-150;
    end
end
track.hdg=90.-track.hdg;
for jj=1:length(track.hdg)
    if track.hdg(jj)>180
        track.hdg(jj)=track.hdg(jj)-360;
    elseif track.hdg(jj)<-180
        track.hdg(jj)=track.hdg(jj)+360;
    end
end
%% Box
% for ii=1:steps
%     if ii<=steps*.25
%         track.stw(ii)=.5;
%         track.hdg(ii)=90;
%     elseif ii<=(steps*.25)+(90*freq)
%         track.stw(ii)=.5;
%         track.hdg(ii)=track.hdg(ii-1)-.1;
%     elseif ii<=steps*.5
%         track.stw(ii)=.5;
%         track.hdg(ii)=0;
%     elseif ii<=(steps*.5)+(90*freq)
%         track.stw(ii)=.5;
%         track.hdg(ii)=track.hdg(ii-1)-.1;
%     elseif ii<=steps*.75
%         track.stw(ii)=.5;
%         track.hdg(ii)=-90;
%     elseif ii<=(steps*.75)+(90*freq)
%         track.stw(ii)=.5;
%         track.hdg(ii)=track.hdg(ii-1)-.1;
%     else 
%         track.stw(ii)=.5;
%         track.hdg(ii)=-180;
%     end
% end
%     track.hdg=90.-track.hdg;
%     for jj=1:length(track.hdg)
%         if track.hdg(jj)>180
%             track.hdg(jj)=track.hdg(jj)-360;
%         elseif track.hdg(jj)<-180
%             track.hdg(jj)=track.hdg(jj)+360;
%         end
%     end
end
