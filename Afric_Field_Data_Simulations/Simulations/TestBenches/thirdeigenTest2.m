%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 50;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Pol_ground = [1;1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);
SNR = 30;
Noise = (10^(-SNR/20))/sqrt(3);
Window_optimal = 10;    %size of window *ones(1,Window_optimal);%

g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));

s1 = alpha*g + beta*v;
s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;

signal_1 = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*(2*pi*rand(3,Window_optimal) - pi));
signal_2 = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*(2*pi*rand(3,Window_optimal) - pi));

%% Fourth Order ESPRIT
h1 = signal_1(1,:);
v1 = signal_1(2,:);
x1 = signal_1(3,:);

h2 = signal_2(1,:);
v2 = signal_2(2,:);
x2 = signal_2(3,:);

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

R1 = S1*S1'/Window_optimal
R11 = pinv(R1)/norm(pinv(R1))
R2 = S1*S2'/Window_optimal;

%         [eigenvec_4,eigenval_4] = eig(pinv(R1_4 + eye_optimal*eye(6))*R2_4);
[eigenvec,eigenval] = eig(R11*R2)
abs(eigenvec)
abs(eigenval)
angle(eigenval)*180/pi