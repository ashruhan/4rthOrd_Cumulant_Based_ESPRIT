%% clear
clear; clc;

%% Initializations

gpol = [1;0];
vpol = [1;1]./sqrt(2);
fig = pi/6;
fiv = pi/3;
window = 100;

g = gpol*randn(1,window);
v = vpol*randn(1,window);
s1 = g + v;
s2 = g*exp(1i*fig) + v*exp(1i*fiv);

%% Second order operations
R11_2 = s1*s1.'/window;
h1_h1 = R11_2(1,1); 
v1_v1 = R11_2(2,2); 
h1_v1 = R11_2(1,2);

R11c_2 = s1*s1'/window;
h1_h1c = R11c_2(1,1); 
v1_v1c = R11c_2(2,2); 
h1_v1c = R11c_2(1,2);

R22_2 = s2*s2.'/window;
h2_h2 = R22_2(1,1); 
v2_v2 = R22_2(2,2); 
h2_v2 = R11_2(1,2);

R12c_2 = s1*s2'/window;
h1_h2c = R12c_2(1,1); 
v1_v2c = R12c_2(2,2); 
h1_v2c = R12c_2(1,2);

%% Fourth Order Operations
R1_4 = (s1.*s1)*(s1.*s1)'/window;
h1_h1_h1c_h1c = R1_4(1,1); 
v1_v1_v1c_v1c = R1_4(2,2); 
h1_v1_h1c_v1c = R1_4(1,2);

R2_4 = (s1.*s1)*(s2.*s2)'/window;
h1_h1_h2c_h2c = R2_4(1,1); 
v1_v1_v2c_v2c = R2_4(2,2); 
h1_v1_h2c_v2c = R2_4(1,2);
%% Cumulant based operaterations
C1(1,1) = h1_h1_h1c_h1c...
        - h1_h1*h1_h1...
        - h1_h1c*h1_h1c...
        - h1_h1c*h1_h1c;

C1(1,2) = h1_v1_h1c_v1c...
        - h1_v1*h1_v1...
        - h1_h1c*v1_v1c...
        - h1_v1c*h1_v1c;

C1(2,1) = C1(1,2);

C1(2,2) = v1_v1_v1c_v1c...
        - v1_v1*v1_v1...
        - v1_v1c*v1_v1c ...
        - v1_v1c*v1_v1c;
%% Cumulant Based Operations
C2(1,1) = h1_h1_h2c_h2c...
        - h1_h1*h2_h2...
        - h1_h1c*h1_h2c...
        - h1_h2c*h1_h2c;
    
C2(1,2) = h1_v1_h2c_v2c...
        - h1_v1*h2_v2...
        - h1_h2c*v1_v2c...
        - h1_v2c*h1_v2c;

C2(2,1) = C2(1,2);

C2(2,2) = v1_v1_v2c_v2c...
        - v1_v1*v2_v2...
        - v1_v2c*v1_v2c...
        - v1_v2c*v1_v2c;
%% ESPRIT 
[eigenvec_4,eigenval_4] = eig(pinv(R1_4)*R2_4)
abs(eigenvec_4)
abs(eigenval_4)
0.5*angle(eigenval_4)*180/pi