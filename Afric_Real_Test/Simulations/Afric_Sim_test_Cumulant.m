%% Afric_Sim_test_Cumulant
clc;clear all;
%% Initializations
g_weight = 0.5; %ground weighting factor
v_weight = 1;   %veg weighting factor
n_weight = 1*10^(-1);

Pol_ground = [1;1;0];  Pol_Cum_ground = [1;0;0]; %ground
Pol_vegitation = [1;0;1]; Pol_Cum_vegitation = [0;1;0]; %vegitation

ground_offset = -pi/3; % ground interferomitry offset
vegitation_offset = -pi/6;    % veg interferomitry offset

Window = 100;    %size of window
Averaged_samples = 50;

est_ground_angle=zeros(Averaged_samples,1); est_vegitation_angle=zeros(Averaged_samples,1);
% mag1=zeros(Averaged_samples,1); mag2=zeros(Averaged_samples,1);


%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Averaged_samples;
    Noise = Averaged_sample*n_weight;
    for i = 1:Window
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        signal_1 = g_weight*g + v_weight*v;
        signal_2 = g_weight*exp(1i*ground_offset)*g + v_weight*exp(1i*vegitation_offset)*v;
        
        AddedNoise_1 = Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        AddedNoise_2 = Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));

        signal_1_Noise = signal_1 + AddedNoise_1;
        signal_2_Noise = signal_2 + AddedNoise_2;
        
        S1 = [signal_1_Noise(1,:).*signal_1_Noise(2,:)
            signal_1_Noise(1,:).*signal_1_Noise(3,:)
            signal_1_Noise(2,:).*signal_1_Noise(3,:)];
        
        S2 = [signal_2_Noise(1,:).*signal_2_Noise(2,:)
            signal_2_Noise(1,:).*signal_2_Noise(3,:)
            signal_2_Noise(2,:).*signal_2_Noise(3,:)];
        
        R1 = S1*S1';
        R2 = S1*S2';
        
        [u,c] = eig(pinv(R1)*R2);
        
        sg = abs(Pol_Cum_ground'*u);
        [~,kk] = sort(sg,'descend');
        est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        %         mag1(Averaged_sample) = mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
        sv = abs(Pol_Cum_vegitation'*u);
        [~,kk] = sort(sv,'descend');
        est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        %         mag2(Averaged_sample) = mag2(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
    end
end
%% Ploting Results

figure(3);title('Ground angle in degrees');
subplot(2,1,1)
plot(est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(-1*ones(Averaged_samples,1)*ground_offset*180/pi,'b.'); %actual Ground Phase
hold off;

% figure(4);
subplot(2,1,2)
title('Vegitation angle in degrees');
plot(est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
hold on;
plot(-1*ones(Averaged_samples,1)*vegitation_offset*180/pi,'r.'); %actual Vegitation phase
hold off;

% figure(2)
% hold on;
% plot(mag1,'b+');
% plot(mag2,'ro');
% hold off;
