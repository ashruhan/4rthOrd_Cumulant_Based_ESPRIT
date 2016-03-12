%% Fuck this right now
clear;clc;
%%
G_O = 20;
V_O = 50;
g_off = G_O*pi/180; % ground interferomitry offset
v_off= V_O*pi/180;    % veg interferomitry offset



Pol_ground = [1;1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);

signal_1 = Pol_ground + Pol_vegitation
signal_2 = Pol_ground*exp(1i*g_off) + Pol_vegitation*exp(1i*v_off)

h1 = signal_1(1,1);
v1 = signal_1(2,1);
x1 = signal_1(3,1);

h2 = signal_2(1,1);
v2 = signal_2(2,1);
x2 = signal_2(3,1);

S1 = [h1.*h1
    v1.*v1
    x1.*x1
    h1.*v1
    h1.*x1
    v1.*x1];

S2 = [h2.*h2
    v2.*v2
    x2.*x2
    h2.*v2
    h2.*x2
    v2.*x2];

R1 = S1*S1';
R2 = S1*S2';

[u,c] = eig(pinv(R1)*R2)