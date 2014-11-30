%% Afric_Sim_test_espirit_NewSort
clc;clear all;
%% Initializations
% Setting up the enviornment
g_weight = 1; %ground weighting factor
v_weight = 1;   %veg weighting factor
n_weight = 1*10^(-3/2);% Added Noise to the System
Pol_ground = [1;1;0]; %Multivatiant Ground
Pol_vegitation = [1;0;1]; %Multivariant Vegitation

Averaging_loop_size = 50;
Signal_samples = 100;    %size of Ensamble Average Window
Noise_samples = 50;
roundz = 0.01;
est_ground_angle = zeros(Noise_samples,1);
est_vegitation_angle = zeros(Noise_samples,1);
est_noise_angle = zeros(Noise_samples,1);
%% Setting offset in degrees then converting to radians
% It help me visualize when entering the degrees instead of radians
ground_offset_degrees = 45; ground_offset_radians = ground_offset_degrees*pi/180; % Ground Interferomitry offset
vegitation_offset_degrees = 60; vegitation_offset_radians = vegitation_offset_degrees*pi/180; % Vegitation Interferomitry offset

%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Noise_samples;
    Noise = Averaged_sample*n_weight;
    for unusedVariable = 1:Averaging_loop_size;
        
        ground =  Pol_ground*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        vegitation =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        
        signal_1 = g_weight*ground + v_weight*vegitation;
        signal_2 = g_weight*exp(1i*ground_offset_radians)*ground + g_weight*exp(1i*vegitation_offset_radians)*vegitation;
        
        AddedNoise_1 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));
        AddedNoise_2 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));
        %% Espirit Algorithm
        
        S1 = signal_1 + AddedNoise_1;
        S2 = signal_2 + AddedNoise_2;
        
        R1 = S1*S1';
        R2 = S1*S2';
        
        A =  pinv(R1)*R2;
        [u,uv] = eig(A);
        
        %% Eigen value Decomposition
        % Used for sorting the eigen values in their respective bins
        
        Eone = mean(abs(u(:,1)*uv(1,1)*pinv(u(:,1)))'); %#ok<*UDIM>
        Etwo = mean(abs(u(:,2)*uv(2,2)*pinv(u(:,2)))');
        Ethree = mean(abs(u(:,3)*uv(3,3)*pinv(u(:,3)))');
        
        [~,kk_one]=sort(Eone,'descend');
        [~,kk_two]=sort(Etwo,'descend');
        [~,kk_three]=sort(Ethree,'descend');
        EthreeNoise = 0; EtwoNoise = 0; EoneNoise = 0;
        %% Determining the Noise signal
        if (abs((Ethree(1,1) - Ethree(1,2)))< roundz)
            if (abs((Ethree(1,1) - Ethree(1,3)))< roundz)
                if (abs((Ethree(1,2) - Ethree(1,3)))< roundz)
                    est_noise_angle(Averaged_sample) = est_noise_angle(Averaged_sample) + angle(uv(3,3))/Averaging_loop_size;
                    EthreeNoise = 1;
                end
            end
        end
        if (abs((Etwo(1,1) - Etwo(1,2)))< roundz)
            if (abs((Etwo(1,1) - Etwo(1,3)))< roundz)
                if (abs((Etwo(1,2) - Etwo(1,3)))< roundz)
                    est_noise_angle(Averaged_sample) = est_noise_angle(Averaged_sample) + angle(uv(2,2))/Averaging_loop_size;
                    EtwoNoise = 1;
                end
            end
        end
        if (abs((Eone(1,1) - Eone(1,2)))< roundz)
            if (abs((Eone(1,1) - Eone(1,3)))< roundz)
                if (abs((Eone(1,2) - Eone(1,3)))< roundz)
                    est_noise_angle(Averaged_sample) = est_noise_angle(Averaged_sample) + angle(uv(1,1))/Averaging_loop_size;
                    EoneNoise = 1;
                end
            end
        end
        %% Determining the Vegitation signal
        
        if ((EthreeNoise == 0)&&(kk_three(1) == 3))
             est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + angle(uv(3,3))/Averaging_loop_size;
        elseif ((EtwoNoise == 0)&&(kk_two(1) == 3))
            est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + angle(uv(2,2))/Averaging_loop_size;
        elseif ((EoneNoise == 0)&&(kk_one(1) == 3))
            est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + angle(uv(1,1))/Averaging_loop_size;
        end
        %% Determining the  Ground signal
        
        
        if ((EthreeNoise == 0)&&(kk_three(1) == 2))
             est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + angle(uv(3,3))/Averaging_loop_size;
        elseif ((EtwoNoise == 0)&&(kk_two(1) == 2))
            est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + angle(uv(2,2))/Averaging_loop_size;
        elseif ((EoneNoise == 0)&&(kk_one(1) == 2))
            est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + angle(uv(1,1))/Averaging_loop_size;
        end
        breakpoint = 1;
    end
end
temp = 1:Noise_samples;
n_snr = 10*log10(1./(temp*n_weight).^2)';
%% Ploting Results
% g_error = -1*est_ground_angle*180/pi - ones(Noise_samples,1)*ground_offset_radians*180/pi;
% v_error = -1*est_vegitation_angle*180/pi - ones(Noise_samples,1)*vegitation_offset_radians*180/pi;
%
% mean(v_error)
% mean(g_error)

% figure(1);title('Ground in blue Veg in red');
% hold on;
% plot(n_snr,g_error,'b')
% plot(n_snr,v_error,'r')
% xlabel('SNR dB');ylabel('error in degrees')
% hold off;

figure(2);title('Ground and Vegitation angle in degrees');
plot(n_snr,-1*est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(n_snr,ones(Noise_samples,1)*ground_offset_radians*180/pi,'b.'); %actual Ground Phase
plot(n_snr,-1*est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(n_snr,ones(Noise_samples,1)*vegitation_offset_radians*180/pi,'r.'); %actual Vegitation phase
plot(n_snr,-1*est_noise_angle*180/pi,'go'); %Estimated Vegitation angle
hold off;


