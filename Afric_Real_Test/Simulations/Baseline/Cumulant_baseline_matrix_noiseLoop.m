%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 1;
n_weight = 1*10^(-2);
g_angle = -pi/16;
v_angle = -pi/24;

delx = 0.5; % Delta x of baseline
baseline = delx; %Physical baseline of two sample arrays

array_size = 6; % the amount of ground samples
Averaged_Signal_samples = 500; %Averaged_Signal_samples in matrix
NoiseITER = 20;

est_ground_angle = zeros(1,NoiseITER);
est_vegitation_angle = zeros(1,NoiseITER);

%% Start of Algorithm
for Noise = 1:NoiseITER;
    n_weight = Noise*n_weight;
    for unusedVariable = 1:Averaged_Signal_samples;
%% Setting up Random Phi that differ across samples
phi_one = 1i*2*pi*rand(1,Averaged_Signal_samples);
phi_two = 1i*2*pi*rand(1,Averaged_Signal_samples);
pos = [1:array_size]'; %#ok<NBRAK>
distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength

%% Element Array size with Signal

X = g_weight.*(exp(-1i*distance*sin(g_angle))*exp(phi_one))... %Ground
    + v_weight.*(exp(-1i*distance*sin(v_angle))*exp(phi_two))...  %Vegetaion
    + n_weight.*sqrt(-2.*log(1-rand(array_size,Averaged_Signal_samples))).*exp(1i*2*pi*rand(array_size,Averaged_Signal_samples));  %Noise

%% Setting up Baselines
% Remember to change the Initial Variable delx if Baseline changes

%S1 = [X(1,:); X(2,:)]; S2 = [X(5,:);X(6,:)]; %Baseline = 4delx
%S1 = [X(2,:);X(3,:)]; S2 = [X(4,:);X(5,:)]; %Baseline = 2delx
S1 = [X(3,:);X(4,:)]; S2 = [X(4,:);X(5,:)]; %Baseline = delx
%S1 = [X(2,:);X(3,:);X(4,:)]; S2 = [X(3,:);X(4,:);X(5,:)]; %Baseline = delx
%S1 = [X(1,:);X(2,:);X(3,:);X(4,:)]; S2 = [X(2,:);X(3,:);X(4,:);X(5,:)]; %Baseline = delx


% S1 = [X(1,:)*X(2,:)';
%     X(1,:)*X(3,:)';
%     X(2,:)*X(3,:)'];
% 
% S2 = [X(4,:)*X(5,:)';
%     X(4,:)*X(6,:)';
%     X(5,:)*X(6,:)']; %Baseline = 3delx
%     
% S1 = [X(1,:)*X(2,:)';
%     X(1,:)*X(3,:)';
%     X(2,:)*X(3,:)'];
% 
% S2 = [X(2,:)*X(3,:)';
%     X(2,:)*X(4,:)';
%     X(3,:)*X(4,:)']; %Baseline = delx

%     S1 = [X(1,1)*X(1,2);
%         X(1,1)*X(1,3);
%         X(1,2)*X(1,3)];
%     
%     S2 = [X(1,3)*X(1,4);
%         X(1,3)*X(1,5);
%         X(1,4)*X(1,5)]; %Baseline = 2delx
%% Espirit Algorithm
R1 = S1*S1';
R2 = S1*S2';
A =  pinv(R1)*R2;

[u,uv] = eig(A);
[~,kk]=sort(angle(diag(uv)),'ascend');

%% Averaging over a matix then averaging that signal in a for loop proportional to the matix rows
% The Smaller Angle ends up haveing The largest angle Eigen Value
est_ground_angle(Noise) = est_ground_angle(Noise) + (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)))/Averaged_Signal_samples; %Estimated ground angle
est_vegitation_angle(Noise) = est_vegitation_angle(Noise) + (asin(angle(uv(kk(2),kk(2)))/(2*pi*baseline)))/Averaged_Signal_samples; %Estimated Vegitation angle
    end
end
%% Ploting Results
figure(1);title('Ground angle in degrees');
plot(est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(ones(1,NoiseITER)*g_angle*180/pi,'b.'); %actual Ground Phase
hold off;

figure(2);title('Vegitation angle in degrees');
plot(est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
hold on;
plot(ones(1,NoiseITER)*v_angle*180/pi,'r.'); %actual Vegitation phase
hold off;