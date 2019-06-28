A = [1 2 3 4;5 6 7 8;9 10 11 12]
B = [13 15 14]
[C,ind] = sort(B, 'ascend')
D = A(ind,:)