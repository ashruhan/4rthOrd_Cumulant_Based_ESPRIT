%% Espirit Algorithm conditional Second order
%% Initializations
clear;clc;
[hh,vv,xx] =  openSARdata();
r = 3; c = 3;
L = r+c+1;
[ylength,xlength] = size(hh.off);

ground_phase_est=zeros(ylength,xlength);
vegitation_phase_est=zeros(ylength,xlength);
ground_mag_est=zeros(ylength,xlength);
vegitation_mag_est=zeros(ylength,xlength);

Pol_ground = [1;1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for col = 1:xlength;
        
        clear('R1','R2','S1','S2','A','uv','kk')
        
        R.r = r;R.row = row;R.ylength = ylength;
        C.c = c;C.col = col;C.xlength = xlength;
        
        if(row<r+1)&&(col<c+1)     % condition 1
            [s1,s2] = Average_Condition1(R,C,hh,vv,xx);
            
        elseif(row<r+1)&&(c+1<col)&&(col<xlength-c)    % condition 2
            [s1,s2] = Average_Condition2(R,C,hh,vv,xx);
            
        elseif(row<r+1)&&(col>xlength-c)   % condition 3
            [s1,s2] = Average_Condition3(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(col<c+1)   % condition 4
            [s1,s2] = Average_Condition4(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(c+1<col)&&(col<xlength-c)  % condition 5
            [s1,s2] = Average_Condition5(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(col>xlength-c)     % condition 6
            [s1,s2] = Average_Condition6(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(col<c+1)    % condition 7
            [s1,s2] = Average_Condition7(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(c+1<col)&&(col<xlength-c)   % condition 8
            [s1,s2] = Average_Condition8(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(col>xlength-c)   % condition 9
            [s1,s2] = Average_Condition9(R,C,hh,vv,xx);
        end
        
        S1(1,:) = s1.h;
        S1(2,:) = s1.v;
        S1(3,:) = s1.x;
        
        S2(1,:) = s2.h;
        S2(2,:) = s2.v;
        S2(3,:) = s2.x;
        
        R1 = S1*S1';
        R2 = S1*S2';
        A = pinv(R1)*R2;
        
        [u,uv] = eig(A);
        
        polCharg = abs(Pol_ground'*u);
        [~,kk] = sort(polCharg,'descend');
        ground_phase_est(row,col) = (uv(kk(1),kk(1)));
%         ground_mag_est(row,col) = abs(uv(kk(1),kk(1)));
%         
%         polCharv = abs(Pol_vegitation'*u);
%         [~,kk] = sort(polCharv,'descend');
%         vegitation_phase_est(row,col) = angle(uv(kk(1),kk(1)));
%         vegitation_mag_est(row,col) = abs(uv(kk(1),kk(1)));
        
    end
end

%% Plotting Results
% figure(1); histg = reshape(g,x*y,1); hist(angle(histg),100);
% figure(2); histv = reshape(v,x*y,1); hist(angle(histv),100)
% figure(3); histn = reshape(n,x*y,1); hist(angle(histn),100);

% figure(4); image(vegitation_phase_est); title('angle(v)');
% figure(5); image(ground_phase_est); title('angle(g)');