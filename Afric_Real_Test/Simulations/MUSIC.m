%% Afric_Sim_test_espirit
clc;clear;
%% Initializations

alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;-1;0]/sqrt(2); 
Pol_vegitation = [1;1;1]/sqrt(3); 

ground_offset = 30*pi/180; % ground interferomitry offset
vegitation_offset =  35*pi/180;    % veg interferomitry offset
samples = 10;
Window = 200;    %size of window
Averaged_samples = 30;

ground_phase_est=zeros(Averaged_samples,1); 
vegitation_phase_est=zeros(Averaged_samples,1);
ground_mag_est=zeros(Averaged_samples,1); 
vegitation_mag_est=zeros(Averaged_samples,1);
SNR = zeros(1,Averaged_samples);
%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques

for Averaged_sample = 1:Averaged_samples;
    
    SNR(Averaged_sample) = Averaged_sample-10;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for n = 1:samples
       
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
                                    
        p1 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        p2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = p1*p1'/Window;
        R2 = p1*p2'/Window;
        
        A = pinv(R1)*R2;
        
        [u,c] = eig(A);
        
        s1 = abs(Pol_ground'*u);
        [~,kk] = sort(s1,'descend');
        ground_phase_est(Averaged_sample) = ground_phase_est(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        ground_mag_est(Averaged_sample) = ground_mag_est(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;

        s2 = abs(Pol_vegitation'*u);
        [~,kk] = sort(s2,'descend');
        vegitation_phase_est(Averaged_sample) = vegitation_phase_est(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        vegitation_mag_est(Averaged_sample) = vegitation_mag_est(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;

    end     
end
%% Plotting Results

% figure(1);
% title('Ground and Vegitation Interferometric Phases');
% xlabel('SNR dB');ylabel('Int Phase (Degrees)');
% hold on;
% plot(SNR,ground_phase_est*180/pi,'ro');
% plot(SNR,-1*ones(Averaged_samples,1)*ground_offset*180/pi,'r+'); %actual Ground Phase
% plot(SNR,vegitation_phase_est*180/pi,'go');
% plot(SNR,-1*ones(Averaged_samples,1)*vegitation_offset*180/pi,'g+'); %actual Vegitation phase
% axis([-10 20 -70 -10]);
% legend('Ground Estimated','Ground Actual','Vegitation Estimated','Vegitaion Actual','Location','northeast')
% hold off

% figure(2)
% title('Ground Vegitation Coherance');
% xlabel('SNR (dB)');ylabel('Maginitude')
% hold on;
% plot(SNR,ground_mag_est,'ro');
% plot(SNR,vegitation_mag_est,'bo');
% legend('Ground Coherance','Vetitaion Coherance','Location','east');
% hold off;