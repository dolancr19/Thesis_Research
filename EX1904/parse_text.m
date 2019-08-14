clear variables
number ='data'+string(86);
text=number+'.txt';
mat=number+'.mat';
[PPS,chan1,chan2,chan3,chan4,index]=textread(text,'%f,%f,%f,%f,%f,%u');
data=[PPS chan1 chan2 chan3 chan4 index];
%[chan2,chan3,index]=textread(text,'%f,%f,%u');
%data=[chan2,chan3,index];
save(mat,'data')
plot(data(:,1:end-1))