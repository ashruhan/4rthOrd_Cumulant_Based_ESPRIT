%% Afric_Sim_test_espirit
%% Initializations
% Setting up the enviornment
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
Pol_ground = [1;1;0]; %Multivatiant Ground
Pol_vegitation = [1;1;1]; %Multivariant Vegitation
ground_offset = -pi/3; % Ground Interferomitry offset
vegitation_offset = pi/6;    % Vegitation Interferomitry offset
Window = 100;    %size of Ensamble Average Window
Averaged_samples = 10;
Noise_weight = 0.1;
Noise = 10; % Added Noise to the System
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
        s2 = exp(1i*ground_offset)*g + exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        
        R1 = s1_Noise*s1_Noise';
        R2 = s1_Noise*s2_Noise';
        
        A=pinv(R1)*R2;
        [u,c] = eig(A);
        
        s1 = abs(Pol_ground'*u);
        [~,kk] = sort(s1,'descend');
        
        phi1(Averaged_sample) = phi1(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        mag1(Averaged_sample) = mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
        s2 = abs(Pol_vegitation'*u);
        [ss,kk] = sort(s2,'descend');
        
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