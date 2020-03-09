clear variables

load calc_new
compare = zeros(10,1000);
for mm = 1:10
    kk = mm + 3*(mm-1);
    ii=1;
        while calc_new(kk,ii) > 0
        ii=ii+1;
        end
    ii = ii-1;
    ind(1,mm) = ii;
    for jj = 1:ii
        % To compare horizontal range error between ray trace and geometric
        %r_uncorrected = calc_new(kk,jj)/tan(calc_new(kk+1,jj));
        %compare(mm,jj) = abs(calc_new(kk+3,jj)-r_uncorrected);
        
        % To compare arrival angle error
        compare(mm,jj) = (180/pi)*abs(calc_new(kk+1,1)-calc_new(kk+1,jj));
        
    end
end
% Plot for position error
%figure
%plot(compare(1,1:ind(1,1)),calc_new(1,1:ind(1,1)),'r',compare(2,1:ind(1,2)),calc_new(5,1:ind(1,2)),'g', compare(3,1:ind(1,3)),calc_new(9,1:ind(1,3)),'b', compare(4,1:ind(1,4)),calc_new(13,1:ind(1,4)),'m', compare(5,1:ind(1,5)),calc_new(17,1:ind(1,5)),'k', compare(6,1:ind(1,6)),calc_new(21,1:ind(1,6)),'r--', compare(7,1:ind(1,7)),calc_new(25,1:ind(1,7)),'g--', compare(8,1:ind(1,8)),calc_new(29,1:ind(1,8)),'b--', compare(9,1:ind(1,9)),calc_new(33,1:ind(1,9)),'m--', compare(10,1:ind(1,10)),calc_new(37,1:ind(1,10)),'k--')
%ylim([0 1600])
%title('Calculated Position Error (Uncorrected) vs Depth with Angles of Transmission between 90 and 45 Degrees')
%set(gca,'YDir','Reverse')
%xlabel('Horizontal Range (m)')
%ylabel('Depth (m)')
%legend('90', '85', '80', '75', '70', '65', '60', '55', '50', '45', 'location', 'northeast')

% Plot for angle error
figure
plot(compare(1,1:ind(1,1)),calc_new(1,1:ind(1,1)),'r',compare(2,1:ind(1,2)),calc_new(5,1:ind(1,2)),'g', compare(3,1:ind(1,3)),calc_new(9,1:ind(1,3)),'b', compare(4,1:ind(1,4)),calc_new(13,1:ind(1,4)),'m', compare(5,1:ind(1,5)),calc_new(17,1:ind(1,5)),'k', compare(6,1:ind(1,6)),calc_new(21,1:ind(1,6)),'r--', compare(7,1:ind(1,7)),calc_new(25,1:ind(1,7)),'g--', compare(8,1:ind(1,8)),calc_new(29,1:ind(1,8)),'b--', compare(9,1:ind(1,9)),calc_new(33,1:ind(1,9)),'m--', compare(10,1:ind(1,10)),calc_new(37,1:ind(1,10)),'k--')
ylim([0 1600])
title('Calculated Arrival Angle Error (Uncorrected) vs Depth with Angles of Transmission between 90 and 45 Degrees')
set(gca,'YDir','Reverse')
xlabel('Angle (degrees)')
ylabel('Depth (m)')
legend('90', '85', '80', '75', '70', '65', '60', '55', '50', '45', 'location', 'northeast')
