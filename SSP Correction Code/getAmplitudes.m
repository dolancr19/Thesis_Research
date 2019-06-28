%% Helper Function: getAmplitudes

function [q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz,czz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% getAmplitudes function: get q, p, A values for ray input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize:
q(i+1)=0;
p(i+1)=0;
A(i+1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INSERT YOUR CODE HERE TO GET p, q, A values   %%%%%%%%%%%%%%%%%

cnn(i+1) = cray(i+1)^2 * czz(i)*xiray(i+1)^2; 
q(i+1) = q(i) + cray(i+1)*sstep*p(i); 
p(i+1) = p(i) - sstep*cnn(i+1)*q(i+1)/cray(i+1)^2; 
 
A(i+1) = 1/(4*pi) *A0 * sqrt(abs((cray(i+1)*cos(theta0))/(rray(i+1)*c0*q(i+1)))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end