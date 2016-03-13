%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);

eye_4 = -0.15;
eye_2 = 0;
G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 1000;
windowoffset = 9;
Window = (windowoffset + 1):100;

ground_angle_4 = zeros(length(Window),1);
ground_angle_2 = zeros(length(Window),1);

vegitation_angle_4 = zeros(length(Window),1);
vegitation_angle_2 = zeros(length(Window),1);

SNR = 10;
Noise = (10^(-SNR/20))/sqrt(3);

%% Matrix Construction
for window = Window; 
    for unusedvariable = 1:Averaged_samples
        
        %% Matrix Construction
        
        g =  Pol_ground*exp(1i*2*pi*rand(1,window));
        v =  Pol_vegitation*exp(1i*2*pi*rand(1,window));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,window))).*exp(1i*2*pi*rand(3,window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,window))).*exp(1i*2*pi*rand(3,window));
        
        %% Second Order Stats
        
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/window;
        R2_2 = S1_2*S2_2'/window;
        
        [eigenvec_2,eigenval_2] =eig((pinv(R1_2 + eye_2*eye(3)))...
            *R2_2,'nobalance');
        
        [~,srt_2]=sort(abs(diag(eigenval_2)),'descend');
        Leig_copol = (eigenval_2(srt_2(1),srt_2(1))^2)....
            *abs(eigenvec_2(1,srt_2(1)))^2....
            + abs(eigenvec_2(2,srt_2(1)))^2;
        
        SLeig_copol = (eigenval_2(srt_2(1),srt_2(1))^2)....
            *abs(eigenvec_2(1,srt_2(2)))^2....
            + abs(eigenvec_2(2,srt_2(2)))^2;
        if (Leig_copol >= SLeig_copol)
            
            ground_angle_2(window - windowoffset) = ground_angle_2(window - windowoffset)....
                + ((ground_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
            
            vegitation_angle_2(window - windowoffset) = vegitation_angle_2(window - windowoffset)....
                + ((vegitation_offset - abs(angle(eigenval_2(srt_2(2),srt_2(2)))))^2)/Averaged_samples;
            
        else
            
            ground_angle_2(window - windowoffset) = ground_angle_2(window - windowoffset)....
                + ((ground_offset - abs(angle(eigenval_2(srt_2(2),srt_2(2)))))^2)/Averaged_samples;
            
            vegitation_angle_2(window - windowoffset) = vegitation_angle_2(window - windowoffset)....
                + ((vegitation_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;          
        end
        
        %% Fourth Order Statistics
        [ Cumulant_11, Cumulant_12 ] = Cumulant( s1_Noise ,s2_Noise,window);
        [eigenvec_4,eigenval_4] = eig((pinv(Cumulant_11+eye_4*eye(6)))...
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
            vegitation_angle_4(window - windowoffset) = vegitation_angle_4(window - windowoffset)....
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
            
            ground_angle_4(window - windowoffset) = ground_angle_4(window - windowoffset)....
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;          
        else
            
            vegitation_angle_4(window - windowoffset) = vegitation_angle_4(window - windowoffset)....
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;
            
            ground_angle_4(window - windowoffset) = ground_angle_4(window - windowoffset)....
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
        end
    end
end

%% Root Mean Square
ground_angle_2 = sqrt(ground_angle_2)*180/pi;
vegitation_angle_2 = sqrt(vegitation_angle_2)*180/pi;

ground_angle_4 = sqrt(ground_angle_4)*180/pi;
vegitation_angle_4 = sqrt(vegitation_angle_4)*180/pi;

%% Plotting Results

figure(1);
title('2nd and 4rth Order ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(Window,(ground_angle_4),'b');
plot(Window,(vegitation_angle_4),'g');
plot(Window,(ground_angle_2),'b');
plot(Window,(vegitation_angle_2),'g');
legend('4rth Order Ground','4rth Order Vegetation','2nd Order Ground','2nd Order Vegetaion','Location','northeast')
hold off
