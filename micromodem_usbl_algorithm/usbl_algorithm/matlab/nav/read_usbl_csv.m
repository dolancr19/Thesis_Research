function [timestamp_dn,az, el,owtt, irphase, t_ms,mfdpow, unitvals] = read_usbl_csv(filename,ver);
% $Id: parse_cst.m 936 2004-04-23 02:29:49Z matt $

if (nargin == 1)
    ver = 1;
end

fd = fopen(filename, 'r');
count = 1;

az = []; % always return something
el = [];
range = [];
timestamp = [];
irphase = [];
t_ms = [];
mfdpow = [];
unitvals = [];
uv = [nan nan nan];mp = nan;tms = nan;irp = [nan nan nan nan];
while 1
    line = fgetl(fd);
    if ~isstr(line), break;end;
    tline =  line(1:findstr(line,'$')-1);
    line(1:findstr(line,'$')-1) = [];
    
    if  (strncmp(line,'$CAMFD', 6))
        tt = line(7:findstr(line,'*')-1);
        tt = strrep(tt, ',',' ');
        [tmp,ok] = str2num(tt);
        if (ok)
            mp = tmp;
        end
    end
    if  (strncmp(line,'$CACST', 6))
        %MFD was associated with packet rather than USBL signal
        mp = nan;
    end
    if  (strncmp(line,'$CAUIR', 6))
        tt = line(7:findstr(line,'*')-1);
        tt = strrep(tt, ',',' ');
        irp = str2num(tt);
    end
    if  (strncmp(line,'$CAUXY', 6))
        tt = line(7:findstr(line,'*')-1);
        tt = strrep(tt, ',',' ');
        [uv,ok] = str2num(tt);
    end
    if  (strncmp(line,'$CAUTM', 6))
        tt = line(7:findstr(line,'*')-1);
        tt = strrep(tt, ',',' ');
        [tms,ok] = str2num(tt);
    end
    if  (strncmp(line,'$CAUSB', 6))
        tt = line(8:findstr(line,'*')-1);
        tt = strrep(tt, ',',' ');
        index = findstr(tt, 'Z');
        tt = [strrep(tt(1:index), '-',' ') tt(index+1:end)];
        tt = strrep(tt, ':',' ');
        tt = strrep(tt, 'T',' ');
        tt = strrep(tt, 'Z',' ');
        [tmp,ok] = str2num(tt);
        
        if (ok)
            yyyy = tmp(1);
            mm = tmp(2);
            dd = tmp(3);
            hr = tmp(4);
            mn = tmp(5);
            ss = tmp(6);
            timestamp_dn(count) = datenum([yyyy, mm,dd,hr,mn,ss]);
            az(count) = tmp(7);
            el(count) = tmp(8);
            owtt(count) = tmp(9);
            mfdpow(count,:) = mp;mp=nan;
            irphase(count,:) = irp; irp = nan.*ones(size(irp));
            t_ms(count) = tms; tms = nan;
            unitvals(count,:) = uv; uv = nan.*ones(length(uv));
            count = count+1;
        end
    end
    
end
