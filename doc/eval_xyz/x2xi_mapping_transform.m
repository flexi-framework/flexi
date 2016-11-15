clear;
P1 = [-1,-1,-1]';
P2 = [+1,-1,-1]';
P3 = [+1,+1,-1]';
P4 = [-1,+1,-1]';
P5 = [-1,-1,+1]';
P6 = [+1,-1,+1]';
P7 = [+1.,+1.,+1.]';
P8 = [-1.,+1,+1]';
m = 1000;
P1 = P1 * m;
P2 = P2 * m;
P3 = P3 * m;
P4 = P4 * m;
P5 = P5 * m;
P6 = P6 * m;
P7 = P7 * m;
P8 = P8 * m;

% transform into unity elem
T(:,1) = 0.5* (P2-P1);
T(:,2) = 0.5* (P4-P1);
T(:,3) = 0.5* (P5-P1);
T_inv = inv(T);

% 
P1T = [-1,-1,-1]';
P2T = [+1,-1,-1]';
P4T = [-1,+1,-1]';
P5T = [-1,-1,+1]';

P3T = T_inv*(P3-P1)-[1,1,1]';
P6T = T_inv*(P6-P1)-[1,1,1]';
P7T = T_inv*(P7-P1)-[1,1,1]';
P8T = T_inv*(P8-P1)-[1,1,1]';

K1 = 1/8*(P1T+P2T+P3T+P4T+P5T+P6T+P7T+P8T);
% for linear case:
%K2 = 1/2* (-P1T+P2T);  
%K3 = 1/2* (-P1T+P4T);
%K4 = 1/2* (-P1T+P5T);
K2 = 1/8*(-P1T+P2T+P3T-P4T-P5T+P6T+P7T-P8T);
K3 = 1/8*(-P1T-P2T+P3T+P4T-P5T-P6T+P7T+P8T);
K4 = 1/8*(-P1T-P2T-P3T-P4T+P5T+P6T+P7T+P8T);

KM(:,1) = K2;
KM(:,2) = K3;
KM(:,3) = K4;
KM_inv = inv(KM);

x = [0.85,0.75,0.65]'*m;
xT = T_inv*(x-P1)-[1,1,1]';
xi = KM_inv*xT-K1;

% Now do newton method:
[F,dF]=NewtonFunc(xi,xT,P1T,P2T,P3T,P4T,P5T,P6T,P7T,P8T);
%for j=1:5
i = 0;
while sum(abs(F)) >= 1E-8
  dF_inv = inv(dF);
  %s=dF\F;
  s=dF_inv*F;
  i=i+1;
  xi=xi-s; % check if s serves as criteria? 
  [F,dF]=NewtonFunc(xi,xT,P1T,P2T,P3T,P4T,P5T,P6T,P7T,P8T);
end
i

% Test
xi

1/8 * (P1*(1-xi(1)) * (1-xi(2)) * (1-xi(3)) ...
         + P2*(1+xi(1)) * (1-xi(2)) * (1-xi(3)) ...
         + P3*(1+xi(1)) * (1+xi(2)) * (1-xi(3)) ...
         + P4*(1-xi(1)) * (1+xi(2)) * (1-xi(3)) ...
         + P5*(1-xi(1)) * (1-xi(2)) * (1+xi(3)) ...
         + P6*(1+xi(1)) * (1-xi(2)) * (1+xi(3)) ...
         + P7*(1+xi(1)) * (1+xi(2)) * (1+xi(3)) ...
         + P8*(1-xi(1)) * (1+xi(2)) * (1+xi(3)))