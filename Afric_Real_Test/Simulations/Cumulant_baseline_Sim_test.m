%% Afric_Sim_test_Cumulant
%% Initializations
Magnitude_1 = 1;
Magnitude_2 = 1/2;

Phase_1 = pi/6;
Phase_2 = pi/4;

array_size = 6;
Window = 100;    %size of window
Averaged_samples = 10;
Noise_weight = 0.1;
X = zeros(1,6);
phi1=zeros(Averaged_samples,1); phi2=zeros(Averaged_samples,1);
mag1=zeros(Averaged_samples,1); mag2=zeros(Averaged_samples,1);


%% Two Element Array size
% A on Case 1
for Averaged_sample = 1:Averaged_samples;
%     Noise = Averaged_sample*Noise_weight;
    for pos = 1:array_size;
        X(pos) =...
              Magnitude_1*exp(1i*2*pi*rand)*exp(1i*Phase_1)...
            + Magnitude_2*exp(1i*2*pi*rand)*exp(1i*Phase_2)...
            + Noise*exp(1i*2*pi*rand);
    end    
        p1=[s1_Noise(1,:).*s1_Noise(2,:)...
            s1_Noise(1,:).*s1_Noise(3,:)...
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        p2=[s2_Noise(1,:).*s2_Noise(2,:)...
            s2_Noise(1,:).*s2_Noise(3,:)...
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1 = p1*p1';
        R2 = p1*p2';
        A=pinv(R1)*R2;
        
        [u,c] = eig(A);
        
        s1=abs(Pol_Cum_ground'*u);
        [~,kk]=sort(s1,'descend');
        phi1(Averaged_sample)=phi1(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        mag1(Averaged_sample)=mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
        s2=abs(Pol_Cum_vegitation'*u);
        [ss,kk]=sort(s2,'descend');
        phi2(Averaged_sample) = phi2(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/Window;
        mag2(Averaged_sample) = mag2(Averaged_sample) + abs(c(kk(1),kk(1)))/Window;
        
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
