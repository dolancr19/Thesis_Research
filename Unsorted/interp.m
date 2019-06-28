%function interp(cmunk,depmat)
load c_munk
load depmat_flat
cGridVal = zeros(length(cmunk),3);
dGridVal = zeros(length(depmat),2);
for ii = 2:length(cmunk)   
    c_prime(ii) = (cmunk(ii,2)-cmunk(ii-1,2))/(cmunk(ii,1)-cmunk(ii-1,1));
    c_2_prime(ii) = (c_prime(ii)-c_prime(ii-1))/(cmunk(ii,1)-cmunk(ii-1,1));
    cGridVal(ii,:) = [cmunk(ii,2), c_prime(ii), c_2_prime(ii)];
end

for jj = 2:length(depmat)
    d_prime(jj) = (depmat(jj,2)-depmat(jj-1,2))/(depmat(jj,1)-depmat(jj-1,1));
    angle(jj) = atand(d_prime(jj));
    dGridVal(jj,:) = [depmat(jj,2), angle(jj)];
end
%end
