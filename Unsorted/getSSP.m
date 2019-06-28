clear all
depvec=(0:7000)';
for ii = 1:length(depvec)
    zbar = 2.*((ii-1)-1300)./1300;
    cvec(ii)=1500.*(1+.00737.*(zbar-1+exp(-zbar)));
    cmunkdeep(ii,:) = [depvec(ii),cvec(ii)];
end