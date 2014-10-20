%% Afric_Sim_test_espirit
%% Initializations
% Setting up the enviornment
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;1;0]./sqrt(2); %Multivatiant Ground
Pol_vegitation = [1;0;1]./sqrt(2); %Multivariant Vegitation
ground_angle = -pi/3; % Ground Interferomitry offset
vegitation_angle = -pi/6;    % Vegitation Interferomitry offset
array_size = 6;
Window = 100;    %size of Ensamble Average Window
Averaged_samples = 20;
Noise_weight = 0.1;% Added Noise to the System
phi1 = zeros(Averaged_samples,1); phi2 = zeros(Averaged_samples,1);
mag1 = zeros(Averaged_samples,1); mag2 = zeros(Averaged_samples,1);

%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Averaged_samples;
    Noise = Averaged_sample*Noise_weight;
    
    phi_one = 1i*2*pi*rand;
    phi_two = 1i*2*pi*rand;
 
    position = 1:array_size;
    distance = 2*pi*position/2;  % element separation of 1/2 wavelength
    
    for i = 1:Window
        
        g_polar =  Pol_ground*exp(phi_one);
        v_polar =  Pol_vegitation*exp(phi_two);
        
        g_polar_angle = g_polar*exp(-1i*distance*sin(ground_angle));
        v_polar_angle = v_polar*exp(-1i*distance*sin(vegitation_angle));
        
        signal = alpha*g_polar_angle + beta*v_polar_angle;
        
        signal_Noise = signal + Noise*sqrt(-2*log(1-rand(3,array_size))).*exp(1i*2*pi*rand(3,array_size));
        
        
        
        
        
        
        R1 = signal_Noise*signal_Noise';
        R2 = signal_Noise*s2_Noise';
        
        A=pinv(R1)*R2;
        [u,c] = eig(A);
        
        sg = abs(Pol_ground'*u);
        [~,kk] = sort(sg,'descend');
        
        phi1(Averaged_sample) = phi1(Averaged_sample) + angle(c(kk(1),kk(1)))/Window;
        mag1(Averaged_sample) = mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
        sv = abs(Pol_vegitation'*u);
        [~,kk] = sort(sv,'descend');
        
        phi2(Averaged_sample) = phi2(Averaged_sample) + angle(c(kk(1),kk(1)))/Window;
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