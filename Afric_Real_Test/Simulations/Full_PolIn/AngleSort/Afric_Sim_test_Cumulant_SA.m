%% Afric_Sim_test_Cumulant_SA
clc;clear all;
%% Initializations
g_weight = 1; %ground weighting factor
v_weight = 1;   %veg weighting factor
n_weight = 1*10^(-3/2);% Added Noise to the System

Pol_ground = [1;1;0]; %ground
Pol_vegitation = [1;0;1];%vegitation

Averaging_loop_size = 50;
Signal_samples = 100;    %size of Ensamble Average Window
Noise_samples = 50;

est_ground_angle=zeros(Noise_samples,1);
est_vegitation_angle=zeros(Noise_samples,1);

%% Setting offset in degrees then converting to radians
% It help me visualize when entering the degrees instead of radians
ground_offset_degrees = 45; ground_offset_radians = ground_offset_degrees*pi/180; % Ground Interferomitry offset
vegitation_offset_degrees = 60; vegitation_offset_radians = vegitation_offset_degrees*pi/180; % Vegitation Interferomitry offset

%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Noise_samples;
    Noise = Averaged_sample*n_weight;
    for unusedVariable = 1:Averaging_loop_size
        
        ground =  Pol_ground*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        vegitation =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        
        signal_1 = g_weight*ground + v_weight*vegitation;
        signal_2 = g_weight*exp(1i*ground_offset_radians)*ground + v_weight*exp(1i*vegitation_offset_radians)*vegitation;
        
        AddedNoise_1 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));
        AddedNoise_2 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));
        
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
        A =  pinv(R1)*R2;
        [u,uv] = eig(A);
        
        temp = 0.5*angle(diag(uv));
        sort_ground_degrees = 0.5*abs(angle(diag(uv)))*180/pi;
        
        [~,kk]=sort(abs(sort_ground_degrees - ground_offset_degrees),'ascend');
        est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + 0.5*angle(uv(kk(1),kk(1)))/Averaging_loop_size;
        
        [~,kk]=sort(abs(sort_ground_degrees - vegitation_offset_degrees),'ascend');
        est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + 0.5*angle(uv(kk(1),kk(1)))/Averaging_loop_size;
    end
end

temp = 1:Noise_samples;
n_snr = 10*log10(1./(temp*n_weight).^2)';
%% Ploting Results
g_error = -1*est_ground_angle*180/pi - ones(Noise_samples,1)*ground_offset_radians*180/pi;
v_error = -1*est_vegitation_angle*180/pi - ones(Noise_samples,1)*vegitation_offset_radians*180/pi;

figure(4)
hold on;
plot(n_snr,g_error,'b')
plot(n_snr,v_error,'r')
hold off;

figure(5);title('Ground and Vegitation angle in degrees');
plot(n_snr,-1*est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(n_snr,ones(Noise_samples,1)*ground_offset_radians*180/pi,'b.'); %actual Ground Phase
plot(n_snr,-1*est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(n_snr,ones(Noise_samples,1)*vegitation_offset_radians*180/pi,'r.'); %actual Vegitation phase
hold off;