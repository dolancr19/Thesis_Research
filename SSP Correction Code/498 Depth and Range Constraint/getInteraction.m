%% Helper Function: getInteraction

function [interact,sray,zray,rray,zetaray,xiray,cray,tau,thetavec,c0,theta0,cz,czz]=getInteraction(r,...
    z,dVal,sray,zray,rray,zetaray,xiray,cray,theta0,cVal,tau,thetavec,c0,cz,czz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% getBottomInteraction function: return interaction %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=length(sray);%current index

%find the depth and bottom angle for your range:
[D,bot_ang]=getDepth(r,dVal);

%initialize vars
interact=false;%no interaction initialized
if (z > D)
    
    %there is a bottom interaction because you want to go lower than
    %the bottom
    %disp('bottom interaction!')
    interact=true;
    %calculate the pathlength needed to reach the bottom exactly
    dsnew=(D-zray(m))/(sin(thetavec(m)));
    dtau=dsnew/cray(m);
    tau(m+1)=tau(m)+dtau;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+1  %%%%%%%%%%%%%%%%
    %%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%
    
    zray(m+1) = zray(m) + dsnew*cray(m)*zetaray(m);
    rray(m+1) = rray(m) + dsnew*cray(m)*xiray(m); 
    
    [cray(m+1), cz(m+1), czz(m+1)] = getCVal(rray(m+1), zray(m+1), cVal);
    
    xiray(m+1) = xiray(m); 
    zetaray(m+1) = zetaray(m) - dsnew/(cray(m)^2)*cz(m); 
    thetavec(m+1) = atan(zetaray(m+1)/xiray(m+1)); 
        
    sray(m+1) = sray(m) + dsnew;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate for next step: going up one step.
    
    % reset c0 and theta0- as though rays are now starting with new angles
    % necessary for when there is a sloped bottom
    c0=cray(m+1);
    theta0= -1*thetavec(m+1)+2*bot_ang;
    dtau=dsnew/cray(m+1);
    tau(m+2)=tau(m+1)+dtau;
    dsnew=(D-zray(m+1))/(sin(thetavec(m+1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+2  %%%%%%%%%%%%%%%%
    %%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%
    
    rray(m+2) = rray(m+1) + dsnew*cray(m+1)*xiray(m+1); 
    zray(m+2) = zray(m+1) + dsnew*cray(m+1)*zetaray(m+1); 
    
    [cray(m+2), cz(m+2), czz(m+2)] = getCVal(rray(m+2), zray(m+2), cVal);
    
    xiray(m+2) = cos(theta0)/c0; 
    zetaray(m+2) = sin(theta0)/c0; 
    thetavec(m+2) = atan(zetaray(m+2)/xiray(m+2)); 
    
    sray(m+2) = sray(m+1) + dsnew;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (z < 0)
    %there is a top interaction
    %disp('top interaction!')
    interact=true;
    
    %calculate the pathlength needed to reach the bottom exactly
    dsnew=abs((zray(m))/(sin(thetavec(m))));
    dtau=dsnew/cray(m);
    tau(m+1)=tau(m)+dtau;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+1  %%%%%%%%%%%%%%%%
    %%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%
    
    zray(m+1) = 0;
    rray(m+1) = rray(m) + dsnew*cray(m)*zetaray(m);
    [cray(m+1), cz(m+1),czz(m+1)] = getCVal(rray(m+1), zray(m+1), cVal);
    zetaray(m+1) = zetaray(m) - dsnew/(cray(m)^2)*cz(m);
    xiray(m+1) = xiray(m); 

    if abs(thetavec(m)) < pi
        thetavec(m+1) = -thetavec(m);
    else 
        thetavec(m+1) = 2*pi - thetavec(m);
    end
    sray(m+1) = sray(m) + dsnew; 
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate for next step: going up one step.
    
    % reset c0 and theta0- as though rays are now starting with new angles
    % necessary for when there is a sloped bottom
    c0=cray(m+1);
    theta0= thetavec(m+1);
    dtau=dsnew/cray(m);
    tau(m+2)=tau(m+1)+dtau;
    dsnew=(0-zray(m+1))/(sin(thetavec(m+1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+2  %%%%%%%%%%%%%%%%
    %%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%
    
    rray(m+2) = rray(m+1);
    zray(m+2) = zray(m+1);
    [cray(m+2), cz(m+2), czz(m+2)] = getCVal(rray(m+2), zray(m+2), cVal);
    
    xiray(m+2) = cos(theta0)/c0;
    zetaray(m+2) = sin(theta0)/c0;
    thetavec(m+2) = atan(zetaray(m+2)/xiray(m+2));
    sray(m+2) = sray(m+1) + dsnew;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
end
