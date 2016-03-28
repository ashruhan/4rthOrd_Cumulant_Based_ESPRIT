%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground_2 = [1;1;0]/sqrt(2);
Pol_ground_4 = [1;1;0;1;0;0]/sqrt(3); %ground
Pol_vegetation_2 = [1;1;1]/sqrt(3);
Pol_vegitation_4 = [1;1;1;1;1;1]/sqrt(6); %vegitation

G_O = 20;
V_O = 60;
ground_offset = G_O*pi/180; % ground interferomitry offset
vegetation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 100;
Window_optimal = 81;    %size of window

eye_distribution = 100;
eye_weight = linspace(-1,1,eye_distribution);

ground_angle_4 = zeros(eye_distribution,1);
ground_angle_2 = zeros(eye_distribution,1);

vegitation_angle_4 = zeros(eye_distribution,1);
vegitation_angle_2 = zeros(eye_distribution,1);

SNR = 10;
Noise = (10^(-SNR/20))/sqrt(3);

%% Matrix Construction
for eye_sample = (1:length(eye_weight));
    for unusedvariable = 1:Averaged_samples
   %% Signal Construction
        
        g =  Pol_ground_2*exp(1i*2*pi*rand(1,Window_optimal));
        v =  Pol_vegetation_2*exp(1i*2*pi*rand(1,Window_optimal));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegetation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        
        index=[1 1
            2 2
            3 3
            1 2
            1 3
            2 3];
        
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        S1_4 = [s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        S2_4 = [s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_4 = S1_4*S1_4'/Window_optimal;
        R2_4 = S1_4*S2_4'/Window_optimal;
        
        R1_2 = S1_2*S1_2'/Window_optimal;
        R2_2 = S1_2*S2_2'/Window_optimal;
 %% Second ORder stats       
        [eigenvec_2,eigenval_2] = eig((pinv(R1_2 + eye_weight(eye_sample)*eye(3)))*R2_2);
        
        polfilter_2 = abs(Pol_ground_2'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        ground_angle_2(eye_sample) = ground_angle_2(eye_sample)...
            + ((ground_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
        
        polfilter_2 = abs(Pol_vegetation_2'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        vegitation_angle_2(eye_sample) = vegitation_angle_2(eye_sample)...
            + ((vegetation_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
        %% Fourth Order Statistics
        [eigenvec_4,eigenval_4] = eig((pinv(R1_4 + eye_weight(eye_sample)*eye(6)))*R2_4);
                
        polfilter_4 = abs(Pol_ground_4'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        ground_angle_4(eye_sample) = ground_angle_4(eye_sample)...
            + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
        
        polfilter_4 = abs(Pol_vegitation_4'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        vegitation_angle_4(eye_sample) = vegitation_angle_4(eye_sample)...
            + ((vegetation_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
        
    end
end

ground_angle_4 = sqrt(ground_angle_4)*180/pi;
ground_angle_2 = sqrt(ground_angle_2)*180/pi;
vegitation_angle_4 = sqrt(vegitation_angle_4)*180/pi;
vegitation_angle_2 = sqrt(vegitation_angle_2)*180/pi;

%% Plotting Results

figure(1);
hold on;
title('4rth Order RMS Error vs Eye Weight');
xlabel('Eye Weight');ylabel('RMS Error (degrees)');
plot(eye_weight,(ground_angle_4),'b');
plot(eye_weight,(vegitation_angle_4),'g');
axis([-1,1,2,15]);
hold off

figure(2);
hold on;
title('2nd Order RMS Error vs Eye Weight');
xlabel('Eye Weight');ylabel('RMS Error (degrees)');
plot(eye_weight,(ground_angle_2),'b');
plot(eye_weight,(vegitation_angle_2),'g');
axis([-1,1,2,15]);
hold off;
