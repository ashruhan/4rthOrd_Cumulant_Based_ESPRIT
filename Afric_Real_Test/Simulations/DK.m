%% Afric_Sim_test_espirit
clc;clear all;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
Pol_ground = [1;-1;0]/sqrt(2); 
% Pol_Cum_ground = [-1;0;0]/2; %ground
Pol_Cum_ground = [1;1;0;-1;0;0]/2; %ground
Pol_vegitation = [1;1;1]/sqrt(3); 
% Pol_Cum_vegitation = [1;1;1]/3; %vegitation
Pol_Cum_vegitation = [1;1;1;1;1;1]/3; %vegitation
ground_offset = pi/6; % ground interferomitry offset
vegitation_offset = pi/5;    % veg interferomitry offset
samples = 10;
Window = 200;    %size of window
Averaged_samples = 30;

phi1=zeros(Averaged_samples,1); phi2=zeros(Averaged_samples,1);
mag1=zeros(Averaged_samples,1); mag2=zeros(Averaged_samples,1);
phi3=zeros(Averaged_samples,1); phi4=zeros(Averaged_samples,1);
mag3=zeros(Averaged_samples,1); mag4=zeros(Averaged_samples,1);

%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
R1=zeros(3,3);
R2=R1;
% R11=zeros(3,3);
% R22=R11;
R11=zeros(6,6);
R22=R11;

% figure(3)
% hold on;

for Averaged_sample = 1:Averaged_samples;
    
    SNR(Averaged_sample)=Averaged_sample-10;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for n = 1:samples
       
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
                        
        
%         p1=[s1_Noise(1,:).*s1_Noise(1,:)
%             s1_Noise(2,:).*s1_Noise(2,:)
%             s1_Noise(3,:).*s1_Noise(3,:)
%             s1_Noise(1,:).*s1_Noise(2,:)
%             s1_Noise(1,:).*s1_Noise(3,:)
%             s1_Noise(2,:).*s1_Noise(3,:)];
%         
%         p2=[s2_Noise(1,:).*s2_Noise(1,:)
%             s2_Noise(2,:).*s2_Noise(2,:)
%             s2_Noise(3,:).*s2_Noise(3,:)
%             s2_Noise(1,:).*s2_Noise(2,:)
%             s2_Noise(1,:).*s2_Noise(3,:)
%             s2_Noise(2,:).*s2_Noise(3,:)];
%         
%         R11 = p1*p1'/Window;
%         R22 = p1*p2'/Window;
                          
        p1=[s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        p2=[s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = p1*p1'/Window;
        R2 = p1*p2'/Window;

%         A=pinv(R11)*R22;
%         
%         [u,c] = eig(A);
%         
%         s1=abs(Pol_Cum_ground'*u);
%         [ss,kk]=sort(s1,'descend');
%         phi1(Averaged_sample)=phi1(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/samples;
%         mag1(Averaged_sample)=mag1(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
%         
% %         plot(SNR(Averaged_sample),0.5*angle(c(kk(1),kk(1)))*180/pi,'b+');
%         
%         s2=abs(Pol_Cum_vegitation'*u);
%         [ss,kk]=sort(s2,'descend');
%         phi2(Averaged_sample) = phi2(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/samples;
%         mag2(Averaged_sample) = mag2(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
        
        A=pinv(R1)*R2;
        
        [u,c] = eig(A);
        
        s1=abs(Pol_ground'*u);
        [ss,kk]=sort(s1,'descend');
        phi3(Averaged_sample)=phi3(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        mag3(Averaged_sample)=mag3(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;

%         plot(SNR(Averaged_sample),angle(c(kk(1),kk(1)))*180/pi,'ro');
        
        s2=abs(Pol_vegitation'*u);
        [ss,kk]=sort(s2,'descend');
        phi4(Averaged_sample) = phi4(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        mag4(Averaged_sample) = mag4(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;

    end     
end
%% Plotting Results

figure(1)
hold on;
% plot(SNR,phi1*180/pi,'b+');
plot(SNR,phi3*180/pi,'ro');
% plot(SNR,phi2*180/pi,'k+');
plot(SNR,phi4*180/pi,'go');
axis([-10 20 -60 0]);
hold off

figure(2)
hold on;
% plot(SNR,mag1,'b+');
plot(SNR,mag3,'ro');
hold off;