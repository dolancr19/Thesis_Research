% script to send usbl pings over the four channel card with built
% in phase shifts
addpath('C:\cygwin\home\ssingh\umodem2_trunk\matlab\nav\jim');
check = 1;
send_pings = 1;
plot_results = 0;
band = 3;
fs = 80000;
probetype = 0; %0=fsk, 1=psk
if (band == 1)
    
    % PSK probe Parameters
    fb_probe = 4000;
    fc_probe = 9760;
    if probetype == 0
        probesyms = 40;
    else
        probesyms = 200;
    end
elseif (band == 2)
    
    % PSK probe Parameters
    fb_probe = 4000;
    fc_probe = 14880;
    if probetype == 0
        probesyms = 40;
    else
        probesyms = 200;
    end
elseif (band == 3)
    
    % PSK probe Parameters
    fb_probe = 4000;
    fc_probe = 25120;
    if probetype == 0
        probesyms = 40;
    else
        probesyms = 200;
    end
elseif (band == 0)
    fb_probe = 5000;
    fc_probe = 25000;
    probesyms = 40;
end

if probetype == 0
    direction = 0;
else
    direction = 1;
end
probe_length = probesyms/fb_probe ;
silence = 0.25;
prepkt = [zeros(floor((silence)*fs)-4000-probe_length*fs,4)];
if direction == 0
    fd1 = fc_probe-fb_probe/2 ;
    fd2 = fc_probe+fb_probe/2 ;
    fb1 = -fb_probe/2;
    fb2 = fb_probe/2;
else
    fd1 = fc_probe+fb_probe/2 ;
    fd2 = fc_probe-fb_probe/2 ;
    fb1 = fb_probe/2;
    fb2 = -fb_probe/2;
end

[y,t]=lfm0(fb1,fb2,probe_length,fb_probe,0);
[ytone1,t]=lfm0(fd1,fd1,0.25,fs,0);
[ytone2,t]=lfm0(fd2,fd2,0.25,fs,0);
[mfdbb,t]=lfm0(fb1,fb2,probe_length,fb_probe*2,0);

c = 1449.2;% sound speed
lambda = c/fc_probe;
vdist = 0.03/2;
if (0.03 > lambda/2)
    sprintf('array spacing is larger than half-wavelength');
end
stepsize = 10;
%azsent = -90:stepsize:90;
azsent = 30;
%elsent = -90:stepsize:90;
elsent = 10;
nreps = 1;
count = 1;
%search regionsize(azerr)
region = [-1 1 -1 1];
% resolution of search
res = 0.02;
% number of iterations
resloops = 4;
%square grid
% array spacing
% equilateral triangle, 2 sensors on x axis
d1 = 0.02286;                    % meters (0.9 in)
array = [ 0 (sqrt(3)/2)*d1 0; d1/2 0 0; -d1/2 0 0];      % x-y pairs
%for beam-to, rotate the array coordinates by 90 degrees
rotM= [cos(pi/2) 0 sin(pi/2);0 1 0;-sin(pi/2) 0 cos(pi/2)];
rotM= [0 0 1;0 1 0;-1 0 0];
v = rotM*array.';
v = v.';
v =  [ 0 0 (sqrt(3)/2)*d1; d1/2 0 0; -d1/2 0 0];      % x-z pairs
%v = [-vdist vdist 0; vdist vdist 0; vdist -vdist 0; -vdist -vdist 0]; ...
%-ones(numsensors,1)*origin;
% diamond grid
%v = [0 vdist 0; vdist 0 0; 0 -vdist 0; -vdist 0 0]; ...
%-ones(numsensors,1)*origin;
ir = [];
az_elm = [];
y_canned = [];
numsensors = size(v,1);
azimsent = [];
elevsent = [];
prepkt = [zeros(floor((silence)*fs)-4000-probe_length*fs,numsensors)];

for az = azsent
    azimuth0_deg = az;
    for el = elsent
        elevation0_deg = el;
        azimuth0_rad = (pi/180)*azimuth0_deg;
        elevation0_rad = (pi/180)*elevation0_deg;
        % normal vector is FROM source TO receiver, hence minus sign
        normal_vector = -[tan(azimuth0_rad) tan(elevation0_rad) 1]';
        unit_normal_vector = normal_vector/sqrt(normal_vector'*normal_vector);
        %origin = [0.015, 0.015, 0];
        %unit_normal_vector = -[-unit_normal_vector(2); unit_normal_vector(3); -unit_normal_vector(1)];
        [iv,phshift(:,count)] = delay_and_phase_shift_bb_signal(v.', fc_probe, unit_normal_vector, y, c, fb_probe, inf);
        [raz(count,:), rel(count,:), xvec(count,:), yvec(count,:), zvec(count,:)] = phaseshift_to_angleofarrival(phshift(:,count), v.', fc_probe, c);
        y_pb = [b2p(iv, [fb_probe,fc_probe]./fs); prepkt];
        
        % check that the phase difference is indeed what was added to the
        % bb signal
        if check
            x = p2b(y_pb,[fb_probe*2,fc_probe]./fs);
            dthresh = 10 ;
            pown = 20 ;
            [ire,ps,dv,p] = mfd1(x,mfdbb',dthresh,pown);
            %[ire,ps,dv,p] = mfd1(y_pb,mf',dthresh,pown);
            [peak_value,peak_sample] = max(abs(ire(:,1)));
            ireval = ire(peak_sample,:);
            %         ireval = cosd(irphase) + sqrt(-1)*sind(irphase);
            ph(:,count) = 180/pi*(phase(ireval)-phase(ireval(1)));
            R = ireval'*ireval;
            ir(:,count) = ireval;
            region2 = bfcoarse(R,v,fc_probe,c,region, res);
            %  [az_elm(count,:),ncalls,i,xline(count),yline(count), zline(count)] = ...
            %	  bfgolden(R,v(:,1:2),fc_probe,c,region2, res);
            [unitvals(count,:),az_elm(count,:), rescnt(count)] = ...
                bf_fast(R,v,fc_probe,c,region2, res, resloops);
            azimsent(count) = az;
            elevsent(count) = el;
        end
        count = count +1;
        
        if (send_pings == 1)
            len = size(y_pb,1);
            y_multi = zeros(len*nreps,numsensors);
            %y_pb(:,[2 3 4]) = zeros(size(y_pb,1),3);
            
            for j=1:nreps
                y_multi((j-1)*len+1:j*len,1) = y_pb(:,1)./max(y_pb(:));
                y_multi((j-1)*len+1:j*len,2) = ytone1.';
                y_multi((j-1)*len+1:j*len,3) = ytone2.';
            end
            
            
            out_multi(real(y_multi)./20,fs);
            %y_canned = [y_canned;y_multi];
            %  pause(1)
        end
    end
end
if check==1
    subplot(2,1,1)
    plot(azsent,az_elm(:,1))
    subplot(2,1,2)
    plot(azsent,az_elm(:,2))
end
return
%plot error maps
azerr = reshape(azimsent-az_elm(:,1).', length(azsent), length(elsent));
elerr = reshape(elevsent-az_elm(:,2).', length(azsent), length(elsent));
%save usblerror_matlabgenerated10loops3el azerr elerr azsent elsent
figure(1);subplot(2,1,1)
pcolor(azsent, elsent, azerr);shading flat; caxis([-10 10]);colorbar
xlabel('azsent')
ylabel('elsent')
title('error in azimuth, matlab generated')
axis('square')
figure(1);subplot(2,1,2)

pcolor(azsent, elsent, elerr);shading flat; caxis([-10 10]);colorbar
xlabel('azsent')
ylabel('elsent')
title('error in elevation, matlab generated')
axis('square')
orient tall
%print('-dpng', 'usblerror_matlabgenerated10loops3el.png')
figure(2)
clf
pcolor(azsent, elsent, sqrt(azerr.^2+elerr.^2));shading flat; caxis([-1 1]);colorbar
xlabel('azsent')
ylabel('elsent')
title('RMS error, matlab generated')
orient portrait
axis('square')
%print('-dpng', 'usblrmserror_matlabgenerated10loops3el.png')
