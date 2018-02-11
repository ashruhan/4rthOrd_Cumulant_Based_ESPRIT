%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 0.1;
v_weight = 1;
n_weight = 0;
g_angle = pi/6;
v_angle = pi/4;
array_size = 6;
Averaged_samples = 100;

%% Two Element Array size
% A on Case 1

phi_one = 1i.*2.*pi.*rand(1,Averaged_samples);
phi_two = 1i.*2*pi.*rand(1,Averaged_samples);
pos = 1:array_size;
distance = 2.*pi.*pos./2;  % element separation of 1/2 wavelength

X = g_weight.*(exp(-1i.*distance.*sin(g_angle))'*exp(phi_one))... %Ground
    + v_weight.*(exp(-1i.*distance.*sin(v_angle))'*exp(phi_two));%...  %Vegetaion
    %+ n_weight.*(sqrt(-2.*log(1-rand(1,Averaged_samples)))*exp(1i*2*pi*rand(1,Averaged_samples)));  %Noise

Y1 = [X(1,:);X(2,:)]; Y2 = [X(5,:);X(6,:)]; %Baseline = 4d

%Y1 = [X(2,:),X(3,:)]; Y2 = [X(4,:),X(5,:)]; %Baseline = 2d
%    Y1 = [X(3,:),X(4,:)]; Y2 = [X(4,:),X(5,:)]; %Baseline = d
%     Y1 = [X(2,:),X(3,:),X(4,:)]; Y2 = [X(3,:),X(4,:),X(5,:)]; %Baseline = d
%     Y1 = [X(1,:);X(2,:);X(3,:);X(4,:)]; Y2 = [X(2,:);X(3,:);X(4,:);X(5,:)]; %Baseline = d

R1 = Y1'*Y1;
R2 = Y1'*Y2;
A = pinv(R1)*R2;

[u,uv] = eig(A);
[~,kk]=sort(angle(diag(uv)),'ascend');
angle(uv(kk(1),kk(1)))*(180/pi)

%% Ploting Results

figure(1)
plot(angle(uv(kk(1),kk(1)))*(180/pi),'b+');

