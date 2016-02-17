%% Second Order ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [ 1; -1 ; 0;]./sqrt(2);
Pol_vegitation = [ 1; 1; 1;]./sqrt(3);

G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 50;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Average_samples = 100;
Window = 50;    %size of window
SNR_samples = 30;

ground_phase = zeros(SNR_samples,1); 
ground_mag = zeros(SNR_samples,1); 

vegitation_phase = zeros(SNR_samples,1);
vegitation_mag = zeros(SNR_samples,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for SNR_sample = 1:SNR_samples;
    
    SNR(SNR_sample)=SNR_sample-10;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Average_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        
        %% ESPRIT Algorithm
        S1 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = S1*S1'/Window;
        R2 = S1*S2'/Window;
        
        [eigenvec,eigenval] = eig(pinv(R1)*R2);
        
        [~,srt]=sort(angle(diag(eigenval)),'descend');
        
        Leig_copol = abs(eigenvec(1,srt(1)))^2 + abs(eigenvec(2,srt(1)))^2;
        SLeig_copol = abs(eigenvec(1,srt(2)))^2 + abs(eigenvec(2,srt(2)))^2;
                
        if (Leig_copol >= SLeig_copol)
            
        ground_phase(SNR_sample) = ground_phase(SNR_sample) + angle(eigenval(srt(1),srt(1)))/Average_samples;
        ground_mag(SNR_sample) = ground_mag(SNR_sample) + abs(eigenval(srt(1),srt(1)))/Average_samples;
        
        vegitation_phase(SNR_sample) = vegitation_phase(SNR_sample) + angle(eigenval(srt(2),srt(2)))/Average_samples;
        vegitation_mag(SNR_sample) = vegitation_mag(SNR_sample) + abs(eigenval(srt(2),srt(2)))/Average_samples;
            
        else
            
        ground_phase(SNR_sample) = ground_phase(SNR_sample) + angle(eigenval(srt(2),srt(2)))/Average_samples;
        ground_mag(SNR_sample) = ground_mag(SNR_sample) + abs(eigenval(srt(2),srt(2)))/Average_samples;
        
        vegitation_phase(SNR_sample) = vegitation_phase(SNR_sample) + angle(eigenval(srt(1),srt(1)))/Average_samples;
        vegitation_mag(SNR_sample) = vegitation_mag(SNR_sample) + abs(eigenval(srt(1),srt(1)))/Average_samples;
            
        end
        
        
%         val1 = abs(Pol_ground'*eigenvec);
%         [~,srt] = sort(val1,'descend');
%         ground_phase(SNR_sample) = ground_phase(SNR_sample) + angle(eigenval(srt(1),srt(1)))/Average_samples;
%         ground_mag(SNR_sample) = ground_mag(SNR_sample) + abs(eigenval(srt(1),srt(1)))/Average_samples;
%         
%         val2 = abs(Pol_vegitation'*eigenvec);
%         [~,srt] = sort(val2,'descend');
%         vegitation_phase(SNR_sample) = vegitation_phase(SNR_sample) + angle(eigenval(srt(1),srt(1)))/Average_samples;
%         vegitation_mag(SNR_sample) = vegitation_mag(SNR_sample) + abs(eigenval(srt(1),srt(1)))/Average_samples;
        
    end
end
%% Plotting Results

figure(1);
title('2nd Order ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(SNR,ground_phase*180/pi,'bo');
plot(SNR,vegitation_phase*180/pi,'go');
legend('2nd Order Ground','2nd Order Vegitaion','Location','northeast')
hold off

figure(2)
title('2nd Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_mag,'bo');
plot(SNR,vegitation_mag,'go');
legend('2nd Order Ground Coherancee','2nd Order Vegitation Coherancee','Location','southeast');
hold off;
