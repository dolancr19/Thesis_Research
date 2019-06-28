clear variables

cable_length = 4000;
output=zeros(10,2);

for ii=1:10
    depth=400*ii;
    theta=asind(depth/cable_length);
    output(ii,1)=theta;
    if ii <= 7
        range = depth/tand(45);
    else
        range = depth/tand(theta);
    end
    output(ii,2)=range;
end

    