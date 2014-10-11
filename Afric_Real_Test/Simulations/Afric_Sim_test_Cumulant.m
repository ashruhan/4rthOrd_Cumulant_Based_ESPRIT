%% Afric_Sim_test_Cumulant
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
Pol_ground = [1;1;0];  Pol_Cum_ground = [1;0;0]; %ground
Pol_vegitation = [1;0;1]; Pol_Cum_vegitation = [0;1;0]; %vegitation
ground_offset = -pi/3; % ground interferomitry offset
vegitation_offset = -pi/6;    % veg interferomitry offset
Window = 100;    %size of window
Averaged_samples = 20;
Noise_weight = 0.1;
phi1=zeros(Averaged_samples,1); phi2=zeros(Averaged_samples,1);
mag1=zeros(Averaged_samples,1); mag2=zeros(Averaged_samples,1);


%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Averaged_samples;
    Noise = Averaged_sample*Noise_weight;
    for i = 1:Window
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        
        p1 = [s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        p2 = [s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1 = p1*p1';
        R2 = p1*p2';
        
        [u,c] = eig(pinv(R1)*R2);
        
        sg = abs(Pol_Cum_ground'*u);
        [~,kk] = sort(sg,'descend');
        phi1(Averaged_sample) = phi1(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        mag1(Averaged_sample) = mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
        sv = abs(Pol_Cum_vegitation'*u);
        [~,kk] = sort(sv,'descend');
        phi2(Averaged_sample) = phi2(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        mag2(Averaged_sample) = mag2(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
    end
end
%% Ploting Results

figure(1)
hold on;
plot(phi1*180/pi,'b+');
plot(phi2*180/pi,'ro');
hold off

figure(2)
hold on;
plot(mag1,'b+');
plot(mag2,'ro');
hold off;
