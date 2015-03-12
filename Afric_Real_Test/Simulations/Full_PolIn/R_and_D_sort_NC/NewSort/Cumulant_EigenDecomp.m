%% Afric_Sim_test_Cumulant_Eigen_Decomp
clc;clear all;
%% Initializations
g_weight = 1; %ground weighting factor
v_weight = 1;   %veg weighting factor
n_weight = 1*10^(-3/2);% Added Noise to the System

Pol_ground = [1;1;0]; %ground
Pol_vegitation = [1;0;1];%vegitation

Averaging_loop_size = 100;
Signal_samples = 100;    %size of Ensamble Average Window
Noise_samples = 50;

est_ground_angle = zeros(Noise_samples,1);
est_vegitation_angle = zeros(Noise_samples,1);
est_noise_angle = zeros(Noise_samples,1);
%% Setting offset in degrees then converting to radians
% It help me visualize when entering the degrees instead of radians
ground_offset_degrees = 50; ground_offset_radians = ground_offset_degrees*pi/180; % Ground Interferomitry offset
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
        EigenVal = 0.5*angle(diag(uv));
        
        %% Eigen value Decomposition
        % Used for sorting the eigen values in their respective bins
        % Also converting to Degrees. This is done because it is easier to
        % use when sorting.
        
        EoneDeg = mean(abs(u(:,1)*uv(1,1)*pinv(u(:,1)))').*180/pi; %#ok<*UDIM>
        EtwoDeg = mean(abs(u(:,2)*uv(2,2)*pinv(u(:,2)))').*180/pi;
        EthreeDeg = mean(abs(u(:,3)*uv(3,3)*pinv(u(:,3)))').*180/pi;
        
        EdprimeOne = diff(EoneDeg,2,2);
        EdprimeTwo = diff(EtwoDeg,2,2);
        EpdrimeThree = diff(EthreeDeg,2,2);
        
        EDecomp_two_diff = [EdprimeOne;
            EdprimeTwo;
            EpdrimeThree];
        
        [~,kk]=sort(EDecomp_two_diff,'ascend');

        est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + EigenVal(kk(3))/Averaging_loop_size;
        est_noise_angle(Averaged_sample) = est_noise_angle(Averaged_sample) + EigenVal(kk(2))/Averaging_loop_size;
        est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + EigenVal(kk(1))/Averaging_loop_size;
    
    end
end
temp = 1:Noise_samples;
n_snr = 10*log10(1./(temp*n_weight).^2)';
%% Ploting Results
g_error = -1*est_ground_angle*180/pi - ones(Noise_samples,1)*ground_offset_radians*180/pi;
v_error = -1*est_vegitation_angle*180/pi - ones(Noise_samples,1)*vegitation_offset_radians*180/pi;

mean(v_error)
mean(g_error)

figure(3);title('Ground in blue Veg in red');
hold on;
plot(n_snr,g_error,'b')
plot(n_snr,v_error,'r')
xlabel('SNR dB');ylabel('error in degrees')
hold off;

figure(4);title('Ground and Vegitation angle in degrees');
plot(n_snr,-1*est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(n_snr,ones(Noise_samples,1)*ground_offset_radians*180/pi,'b.'); %actual Ground Phase
plot(n_snr,-1*est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(n_snr,ones(Noise_samples,1)*vegitation_offset_radians*180/pi,'r.'); %actual Vegitation phase
% plot(n_snr,-1*est_noise_angle*180/pi,'go'); %Estimated Vegitation angle
hold off;
