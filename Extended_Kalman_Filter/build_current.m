function [current]=build_current(steps, set0, drift0, variation)
current.set=zeros(steps,1);
current.drift=zeros(steps,1);

%% Constant current
if variation==0
    current.set(:)=set0;
    current.drift(:)=drift0;
else

%% Variable current
    for ii=1:steps
        current.set(ii)=variation*sin(ii*1.454e-4)+set0;
    end
    current.drift(:)=drift0;
end
end
