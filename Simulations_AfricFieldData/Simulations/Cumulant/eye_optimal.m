%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground_2 = [1;1;0]/sqrt(2);

Pol_vegitation_2 = [1;1;1]/sqrt(3);

G_O = 30;
V_O = 50;
ground_offset = G_O*pi/180; % ground interferomitry offset
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 1000;
Window_optimal = 81;    %size of window

eye_distribution = 100;
eye_weight = linspace(-2,2,eye_distribution);

ground_angle_4 = zeros(eye_distribution,1);
ground_angle_2 = zeros(eye_distribution,1);

vegitation_angle_4 = zeros(eye_distribution,1);
vegitation_angle_2 = zeros(eye_distribution,1);

SNR = 10;
Noise = (10^(-SNR/20))/sqrt(3);

%% Matrix Construction
for eye_sample = (1:length(eye_weight));   
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground_2*exp(1i*2*pi*rand(1,Window_optimal));
        v =  Pol_vegitation_2*exp(1i*2*pi*rand(1,Window_optimal));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        
        %% Fourth Order Statistics
        [ Cumulant_11, Cumulant_12, Cumulant_22 ] = Cumulant( s1_Noise ,s2_Noise,Window_optimal );
        [eigenvec_4,eigenval_4] = eig((pinv(Cumulant_11 + eye_weight(eye_sample)*eye(6)))...
            *Cumulant_12,'nobalance');

        [~,srt_4] = sort(sqrt(abs(diag(eigenval_4))),'descend');
        
        
        LeigTemp  = (abs(eigenval_4(srt_4(1),srt_4(1))))^2....
            *(abs(eigenvec_4(3,srt_4(1)))^2....
            + abs(eigenvec_4(5,srt_4(1)))^2....
            + abs(eigenvec_4(6,srt_4(1)))^2);
        
        SLeigTemp = (abs(eigenval_4(srt_4(2),srt_4(2))))^2....
            *(abs(eigenvec_4(3,srt_4(2)))^2....
            + abs(eigenvec_4(5,srt_4(2)))^2....
            + abs(eigenvec_4(6,srt_4(2)))^2);
        
        if LeigTemp >= SLeigTemp
            vegitation_angle_4(eye_sample) = vegitation_angle_4(eye_sample)....
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
            
            ground_angle_4(eye_sample) = ground_angle_4(eye_sample)....
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;
        else
            
            vegitation_angle_4(eye_sample) = vegitation_angle_4(eye_sample)....
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;
            
            ground_angle_4(eye_sample) = ground_angle_4(eye_sample)....
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;        
        end 
        %% Second Order ESPRIT
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window_optimal;
        R2_2 = S1_2*S2_2'/Window_optimal;
        
        [eigenvec_2,eigenval_2] =eig((pinv(R1_2 + eye_weight(eye_sample)*eye(3)))...
            *R2_2,'nobalance');
        
        [~,srt_2]=sort(abs(diag(eigenval_2)),'descend');
        Leig_copol = (eigenval_2(srt_2(1),srt_2(1))^2)....
            *abs(eigenvec_2(1,srt_2(1)))^2....
            + abs(eigenvec_2(2,srt_2(1)))^2;
        
        SLeig_copol = (eigenval_2(srt_2(1),srt_2(1))^2)....
            *abs(eigenvec_2(1,srt_2(2)))^2....
            + abs(eigenvec_2(2,srt_2(2)))^2;
        if (Leig_copol >= SLeig_copol)         
            ground_angle_2(eye_sample) = ground_angle_2(eye_sample)....
                + ((ground_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
            
            vegitation_angle_2(eye_sample) = vegitation_angle_2(eye_sample)....
                + ((vegitation_offset - abs(angle(eigenval_2(srt_2(2),srt_2(2)))))^2)/Averaged_samples;
            
        else       
            ground_angle_2(eye_sample) = ground_angle_2(eye_sample)....
                + ((ground_offset - abs(angle(eigenval_2(srt_2(2),srt_2(2)))))^2)/Averaged_samples;
            
            vegitation_angle_2(eye_sample) = vegitation_angle_2(eye_sample)....
                + ((vegitation_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
        end
        
    end
end

ground_angle_4 = sqrt(ground_angle_4)*180/pi;
ground_angle_2 = sqrt(ground_angle_2)*180/pi;
vegitation_angle_4 = sqrt(vegitation_angle_4)*180/pi;
vegitation_angle_2 = sqrt(vegitation_angle_2)*180/pi;

%% Plotting Results

figure(1);
hold on;
title('RMS Error 4rth Order')
plot(eye_weight,(ground_angle_4),'b');
plot(eye_weight,(vegitation_angle_4),'g');
hold off

figure(2);
hold on;
title('RMS Error 2nd Order')
plot(eye_weight,(ground_angle_2),'b');
plot(eye_weight,(vegitation_angle_2),'g');
hold off;
