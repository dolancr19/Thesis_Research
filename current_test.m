clear variables
clc

t=43200;
current=zeros(1,t);
for ii=1:t
    current(ii)=.1*sin(ii*1.454e-4)+.5;
end

figure
plot(1:t,current)