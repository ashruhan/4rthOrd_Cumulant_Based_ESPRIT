index=[1 1
    2 2
    3 3
    1 2
    1 3
    2 3];

S1_2 = [s1_Noise(1,:)
    s1_Noise(2,:)
    s1_Noise(3,:)];

S2_2 = [s2_Noise(1,:)
    s2_Noise(2,:)
    s2_Noise(3,:)];

S1_4 = [s1_Noise(1,:).*s1_Noise(1,:)
    s1_Noise(2,:).*s1_Noise(2,:)
    s1_Noise(3,:).*s1_Noise(3,:)
    s1_Noise(1,:).*s1_Noise(2,:)
    s1_Noise(1,:).*s1_Noise(3,:)
    s1_Noise(2,:).*s1_Noise(3,:)];

S2_4 = [s2_Noise(1,:).*s2_Noise(1,:)
    s2_Noise(2,:).*s2_Noise(2,:)
    s2_Noise(3,:).*s2_Noise(3,:)
    s2_Noise(1,:).*s2_Noise(2,:)
    s2_Noise(1,:).*s2_Noise(3,:)
    s2_Noise(2,:).*s2_Noise(3,:)];

R11 = S1_4*S1_4'/Window_optimal;
R22 = S1_4*S2_4'/Window_optimal;

R1 = S1_2*S1_2'/Window_optimal;
R2 = S1_2*S2_2'/Window_optimal;

R10 = S1_2*S1_2.'/Window_optimal;
R20 = S2_2*S2_2.'/Window_optimal;

for m = 1:6 
    for n = 1:6
        
        R11(m,n)= R11(m,n) - R10(index(m,1),index(m,2))*R10(index(n,1),index(n,2))...
            - R1(index(m,1),index(n,1))*R1(index(m,2),index(n,2))...
            - R1(index(m,1),index(n,2))*R1(index(m,2),index(n,1));
        
        R22(m,n)= R22(m,n) - R10(index(m,1),index(m,2))*R20(index(n,1),index(n,2))...
            - R2(index(m,1),index(n,1))*R2(index(m,2),index(n,2))...
            - R2(index(m,1),index(n,2))*R2(index(m,2),index(n,1));
        
    end 
end