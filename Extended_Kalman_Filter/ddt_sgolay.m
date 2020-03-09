function dpxdtp = ddt_sgolay(x,order,window,dt,p)
% dpxdtp = ddt_sgolay(x,order,window,dt,p)
% Wrapper aroung sgolay for computing pth derivative of a noisy signal.
% Applies only to regularly spaced data.
%
% Revision History
% 2017-12-20    mvj    Created.


[~,g] = sgolay(order,window);
dpxdtp = conv(x,factorial(p)/(-dt)^p * g(:,p+1), 'same');
