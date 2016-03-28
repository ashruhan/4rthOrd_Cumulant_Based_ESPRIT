%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;-1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);

ground_offset = 30*pi/180; % ground interferomitry offset
vegitation_offset = 50*pi/180;    % veg interferomitry offset
%
Averaged_samples = 1000;

Window = 81;    %size of window

eye_samples = 100;
%
ground_phase_4=zeros(eye_samples,1);
ground_mag_4=zeros(eye_samples,1);

vegitation_phase_4=zeros(eye_samples,1);
vegitation_mag_4=zeros(eye_samples,1);

forth_order = 6;
sort_ground_4 = zeros(forth_order,1);

second_order = 3;
sort_ground_2 = zeros(second_order,1);
sort_vegetation_2 = zeros(second_order,1);
ground_phase_est_fourth = zeros(1,eye_samples);
% SNR = zeros(1,SNR_samples);
SNR = 10;
eye_w = linspace(0,1,eye_samples);
%% Matrix Construction
for eye_sample_sample = (1:length(eye_w));
    
    Noise = (10^(-SNR/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(second_order,Window))).*exp(1i*2*pi*rand(second_order,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(second_order,Window))).*exp(1i*2*pi*rand(second_order,Window));
        
        
        %% Fourth Order ESPRIT
        S1_4=[s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        S2_4=[s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_4 = S1_4*S1_4'/Window;
        R2_4 = S1_4*S2_4'/Window;
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4 + eye_w(eye_sample_sample)*eye(6))*R2_4);
        
        for i = 1:forth_order
            sort_ground_4(i) = (abs(eigenval_4(i,i)))*(abs(eigenvec_4(1,i))^2 + abs(eigenvec_4(2,i))^2 + abs(eigenvec_4(4,i))^2);
        end       
        
        [~,srt_g_4] = sort(sort_ground_4,'descend');
        
        ground_phase_est_fourth(eye_sample_sample) = ground_phase_est_fourth(eye_sample_sample)...
            + ((ground_offset - abs(0.5*angle(eigenval_4(srt_g_4(1),srt_g_4(1)))))^2)/Averaged_samples;
        
    end
end
ground_phase_est_fourth_mnsqrd = sqrt(ground_phase_est_fourth)*180/pi;
%% Plotting Results

figure(2);title('RMS Error (dB)')
plot((ground_phase_est_fourth_mnsqrd),'black');
xlabel('SNR');ylabel('RMS Error(dB)')
