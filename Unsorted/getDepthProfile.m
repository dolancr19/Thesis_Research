clear all
D = 7000;
depvec=(0:D)';
for ii = 1:length(depvec)
    depmatdeep(ii,:) = [(ii-1),D];
end