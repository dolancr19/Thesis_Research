function [xx,yy]=northriver_xy(lat,lon)
%returns the along-river coordinate (x,y) of a site (lat,lon) in north river
%x=0 at the mouth and decreases upstream (x<0)
%y=0 on main stream

% x shift
mouth_pos=[-70.71158 42.16191];

load north_river_xy.mat;

[~,ii]=min(abs( (lonb-mouth_pos(1))+1i*(latb-mouth_pos(2)) ));

lat=lat(:); lon=lon(:);

xx=nan*ones(size(lat)); yy=xx; arr=xx; beta=xx;
ind=ones(size(lat));

dis_b=dis_b-dis_b(end); %% dis=0 at the mouth and then decreases upstream

xshift=dis_b(ii);
dis_b=dis_b-xshift;

[x,y]=latlon2utm(lat,lon);

x=x-x0; y=y-y0;

%%
xaxiss=diff(xb)+1i*diff(yb);

betas=angle(xaxiss/(1+1i*0))/pi*180; %angle from old xaxis to new xaxis
    %i.e. from Cartesian coordinate to local (along-river) coordinate
    
%% first round
lenx=length(x);
lenxb=length(xb);

for k=1:lenx
    dis_arr=(x(k)-xb)+1i*(y(k)-yb);
    [~,ind(k)]=min( abs( dis_arr ) );
    
    if ind(k)>=lenxb
        ind(k)=lenxb-1;
    end
    arr(k)=dis_arr(ind(k));
end
                 
xaxis1=xaxiss(ind);

alpha=angle(xaxis1./arr)/pi*180;  %% angle from arr to xaxis1 (new xaxis)

ind( abs(alpha)>90 )=ind( abs(alpha)>90 )-1;    %%% 
ind=max(ind,1);

%% second round
for k=1:lenx
    dis_arr=(x(k)-xb)+1i*(y(k)-yb);
    arr(k)=dis_arr(ind(k));
end
% xaxis1=xaxiss(ind);
% alpha=angle(xaxis1./arr)/pi*180;
beta=betas(ind);

%%
[xx,yy]=rot(real(arr),imag(arr),-beta);
xx=xx+dis_b(ind);

return


function [xr,yr]=rot(x,y,theta,ya)
if (nargin==2)||(length(x)~=length(y))
    if nargin==3;
        ya=theta;
        narg=4;
    end
    theta=y;
    y=imag(x);
    x=real(x);
    comp=1;
    narg=2;
else
    comp=0;
    narg=nargin;
end

if narg==4
    theta=-atan2(ya,theta)*180/pi;
end
if imag(theta)~=0
    theta=-atan2(imag(theta),real(theta))*180/pi;
end
costheta=cos(theta/180*pi);
sintheta=sin(theta/180*pi);

if length(theta)==1
    xr=x*costheta-y*sintheta;
    yr=x*sintheta+y*costheta;
else
    xr=x.*costheta-y.*sintheta;
    yr=x.*sintheta+y.*costheta;
end

if comp==1;
    xr=xr+1i*yr;
    yr=[];
end
return;

