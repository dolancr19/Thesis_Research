for jj=1:30
    kk=((jj-1)*20)+1;
    A=ones(2,20);
    A(2,:)=data(7,kk:(jj*20));
A=A';

b=data(8,kk:(jj*20));
b=b';

x_plus=pinv(A)*b;
new_y=A*x_plus;
smooth_y(kk:(jj*20))=new_y;
end

plot(data(7,1:600),smooth_y(1,:));%data(7,1:600),data(8,1:600))