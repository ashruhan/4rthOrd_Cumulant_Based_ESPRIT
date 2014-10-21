%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 0.1;
n_weight = 0;
g_angle = pi/24;
v_angle = pi/16;
array_size = 6;
Averaged_samples = 100;
delx=0.5;

%% Setting up Random Phi that differ across samples

phi_one = 1i*2*pi*rand(1,Averaged_samples);
phi_two = 1i*2*pi*rand(1,Averaged_samples);
pos = 1:array_size;
pos = pos';
distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength

%% Two Element Array size
% A on Case 1

X = g_weight.*(exp(-1i*distance*sin(g_angle))*exp(phi_one))... %Ground
    + v_weight.*(exp(-1i*distance*sin(v_angle))*exp(phi_two))...  %Vegetaion
    + n_weight.*sqrt(-2.*log(1-rand(array_size,Averaged_samples))).*exp(1i*2*pi*rand(array_size,Averaged_samples));  %Noise


Y1 = [X(1,:); X(2,:)]; Y2 = [X(5,:);X(6,:)]; %Baseline = 4d
baseline=4*delx;

% Y1 = [X(1,2),X(1,3)]; Y2 = [X(1,4),X(1,5)]; %Baseline = 2d
%    Y1 = [X(1,3),X(1,4)]; Y2 = [X(1,4),X(1,5)]; %Baseline = d
%     Y1 = [X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,3),X(1,4),X(1,5)]; %Baseline = d
%     Y1 = [X(1,1),X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,2),X(1,3),X(1,4),X(1,5)]; %Baseline = d

R1 = Y1*Y1';
R2 = Y1*Y2';
A =  pinv(R1)*R2;

[u,uv] = eig(A);
[~,kk]=sort(angle(diag(uv)),'ascend');

phase_one = asin(angle(uv(kk(1),kk(1)))/2/pi/baseline);
phase_two = asin(angle(uv(kk(2),kk(2)))/2/pi/baseline);
%% Ploting Results
figure(1)
plot(phase_one*180/pi,'b+')
hold on;
plot(phase_two*180/pi,'r+')
hold off;

