%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations

Pol_ground = [1;0]./sqrt(1);
Pol_vegitation = [1;-1]./sqrt(2);

G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Window_optimal = 1000;

% g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
% v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));


g =  Pol_ground*(exp(1i*2*pi*rand(1,Window_optimal)));
v =  Pol_vegitation*(exp(1i*2*pi*rand(1,Window_optimal)));


s1 = g + v;
s2 = exp(1i*ground_offset)*g + exp(1i*vegitation_offset)*v;

%% Forming the Six Arrays
S1_2 = [s1(1,:)%h
    s1(2,:)];%v

S2_2 = [s2(1,:)%h
    s2(2,:)];%v

R11_2 = S1_2*S1_2.'/Window_optimal;
E_h1_h1 = R11_2(1,1); 
E_h1_v1 = R11_2(1,2);
E_v1_v1 = R11_2(2,2);

R11c_2 = S1_2*S1_2'/Window_optimal;
E_h1_h1c = R11c_2(1,1); 
E_h1_v1c = R11c_2(1,2);
E_v1_h1c = R11c_2(2,1); 
E_v1_v1c = R11c_2(2,2);

R22_2 = S2_2*S2_2.'/Window_optimal;
E_h2_h2 = R22_2(1,1); 
E_h2_v2 = R22_2(1,2);
E_v2_v2 = R22_2(2,2);

R12c_2 = S1_2*S2_2'/Window_optimal;
E_h1_h2c = R12c_2(1,1); 
E_h1_v2c = R12c_2(1,2);
E_v1_h2c = R12c_2(2,1); 
E_v1_v2c = R12c_2(2,2);

S1_4 = [s1(1,:).*s1(1,:)%hh
    s1(2,:).*s1(2,:)%vv
    s1(1,:).*s1(2,:)];%hv

S2_4 = [s2(1,:).*s2(1,:)%hh
    s2(2,:).*s2(2,:)%vv
    s2(1,:).*s2(2,:)];%hv

R1_4 = S1_4*S1_4'/Window_optimal;
R2_4 = S1_4*S2_4'/Window_optimal;

%% Forming Cumulant One Matrix
Cumulant1_4(1,1) = R1_4(1,1) - E_h1_h1*conj(E_h1_h1)...
                             - E_h1_h1c*E_h1_h1c...
                             - E_h1_h1c*E_h1_h1c;
                        
Cumulant1_4(1,2) = R1_4(1,2) - E_h1_h1*conj(E_v1_v1)...
                             - E_h1_v1c*E_h1_v1c...
                             - E_h1_v1c*E_h1_v1c;

Cumulant1_4(1,3) = R1_4(1,3) - E_h1_h1*conj(E_h1_v1)...
                             - E_h1_h1c*E_h1_v1c...
                             - E_h1_v1c*E_h1_h1c;
%%
Cumulant1_4(2,1) = R1_4(2,1) - E_v1_v1*conj(E_h1_h1)...
                             - E_v1_h1c*E_v1_h1c...
                             - E_v1_h1c*E_v1_h1c;

Cumulant1_4(2,2) = R1_4(2,2) - E_v1_v1*conj(E_v1_v1)...
                             - E_v1_v1c*E_v1_v1c...
                             - E_v1_v1c*E_v1_v1c;
                        
Cumulant1_4(2,3) = R1_4(2,3) - E_v1_v1*conj(E_h1_v1)...
                             - E_v1_h1c*E_v1_v1c...
                             - E_v1_v1c*E_v1_h1c;
%%
Cumulant1_4(3,1) = R1_4(3,1) - E_h1_v1*conj(E_h1_h1)...
                             - E_h1_h1c*E_v1_h1c...
                             - E_h1_h1c*E_v1_h1c;
                        
Cumulant1_4(3,2) = R1_4(3,2) - E_h1_v1*conj(E_v1_v1)...
                             - E_h1_v1c*E_v1_v1c...
                             - E_h1_v1c*E_v1_v1c;
                        
Cumulant1_4(3,3) = R1_4(3,3) - E_h1_v1*conj(E_h1_v1)...
                             - E_h1_h1c*E_v1_v1c...
                             - E_h1_v1c*E_v1_h1c;

%% Forming Cumulant Two Matrix
Cumulant2_4(1,1) = R2_4(1,1) - E_h1_h1*conj(E_h2_h2)...
                             - E_h1_h2c*E_h1_h2c...
                             - E_h1_h2c*E_h1_h2c;

Cumulant2_4(1,2) = R2_4(1,2) - E_h1_h1*conj(E_v2_v2)...
                             - E_h1_v2c*E_h1_v2c...
                             - E_h1_v2c*E_h1_v2c;
                        
Cumulant2_4(1,3) = R2_4(1,3) - E_h1_h1*conj(E_h2_v2)...
                             - E_h1_h2c*E_h1_v2c...
                             - E_h1_v2c*E_h1_h2c;
%%
Cumulant2_4(2,1) = R2_4(2,1) - E_v1_v1*conj(E_h2_h2)...
                             - E_v1_h2c*E_v1_h2c...
                             - E_v1_h2c*E_v1_h2c;
                        
Cumulant2_4(2,2) = R2_4(2,2) - E_v1_v1*conj(E_v2_v2)...
                             - E_v1_v2c*E_v1_v2c...
                             - E_v1_v2c*E_v1_v2c;

Cumulant2_4(2,3) = R2_4(2,3) - E_v1_v1*conj(E_h2_v2)...
                             - E_v1_h2c*E_v1_v2c...
                             - E_v1_v2c*E_v1_h2c;
%%
Cumulant2_4(3,1) = R2_4(3,1) - E_h1_v1*conj(E_h2_h2)...
                             - E_h1_h2c*E_v1_h2c...
                             - E_h1_h2c*E_v1_h2c;
                        
Cumulant2_4(3,2) = R2_4(3,2) - E_h1_v1*conj(E_v2_v2)...
                             - E_h1_v2c*E_v1_v2c...
                             - E_h1_v2c*E_v1_v2c;
                        
Cumulant2_4(3,3) = R2_4(3,3) - E_h1_v1*conj(E_h2_v2)...
                             - E_h1_h2c*E_v1_v2c...
                             - E_h1_v2c*E_v1_h2c;

%% ESPRIT Algorithm
[eigenvec_4,eigenval_4] = eig(pinv(Cumulant1_4)*Cumulant2_4);

abs(eigenvec_4)
abs(eigenval_4)
0.5*angle(eigenval_4)*180/pi
