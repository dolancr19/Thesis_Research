clear variables
load plotstuff

figure
plot(plotstuff(:,1),plotstuff(:,2),'r', plotstuff(1:622,3), plotstuff(1:622,4),'b', plotstuff(1:602,5), plotstuff(1:602,6),'g')
title('Calculated Position Error vs Depth with 45, 75 and 90 Degree Angles of Transmission')
xlabel('Horizontal Error (m)')
ylabel('Depth (m)')
legend('45', '75', '90', 'location', 'southeast')