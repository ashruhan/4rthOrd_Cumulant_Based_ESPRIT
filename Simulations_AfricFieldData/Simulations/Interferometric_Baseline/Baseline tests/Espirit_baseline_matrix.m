%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 1;
n_weight = 0;
g_angle = pi/24;
v_angle = pi/16;
array_size = 6;
Averaged_samples = 100;
delx=0.5;
baseline=2*delx;
%% Setting up Random Phi that differ across samples

phi_one = 1i*2*pi*rand(1,Averaged_samples);
phi_two = 1i*2*pi*rand(1,Averaged_samples);
pos = [1:array_size]';
distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength

%% Two Element Array size
% A on Case 1

X = g_weight.*(exp(-1i*distance*sin(g_angle))*exp(phi_one))... %Ground
    + v_weight.*(exp(-1i*distance*sin(v_angle))*exp(phi_two))...  %Vegetaion
    + n_weight.*sqrt(-2.*log(1-rand(array_size,Averaged_samples))).*exp(1i*2*pi*rand(array_size,Averaged_samples));  %Noise


%S1 = [X(1,:); X(2,:)]; S2 = [X(5,:);X(6,:)]; %Baseline = 4d
S1 = [X(2,:);X(3,:)]; S2 = [X(4,:);X(5,:)]; %Baseline = 2d
%S1 = [X(3,:),X(4,:)]; S2 = [X(4,:),X(5,:)]; %Baseline = d
%S1 = [X(2,:),X(3,:),X(4,:)]; S2 = [X(3,:),X(4,:),X(5,:)]; %Baseline = d
%S1 = [X(1,:),X(2,:),X(3,:),X(4,:)]; S2 = [X(2,:),X(3,:),X(4,:),X(5,:)]; %Baseline = d

R1 = S1*S1';
R2 = S1*S2';
A =  pinv(R1)*R2;

[u,uv] = eig(A);
[~,kk]=sort(angle(diag(uv)),'ascend');

est_ground_angle = asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)); 
est_vegitation_angle = asin(angle(uv(kk(2),kk(2)))/(2*pi*baseline));
%% Ploting Results
figure(1)
plot(est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(g_angle*180/pi,'b.'); %actual Ground Phase
plot(v_angle*180/pi,'r.'); %actual Vegitation phase
hold off;

