%% Espirit Algorithm Second order
clear;clc;
%% Initializations
fid(1)=fopen('ref_hh.dat','r');fid(2)=fopen('ref_vv.dat','r');
fid(3)=fopen('ref_xx.dat','r');fid(4)=fopen('offset_hh.dat','r');
fid(5)=fopen('offset_vv.dat','r');fid(6)=fopen('offset_xx.dat','r');

h1=fread(fid(1),[730 1662],'single');h2=fread(fid(4),[730 1662],'single');
v1=fread(fid(2),[730 1662],'single');v2=fread(fid(5),[730 1662],'single');
x1=fread(fid(3),[730 1662],'single');x2=fread(fid(6),[730 1662],'single');

for close = 1:6;
    fclose(fid(close));
end

hh.ref=h1(1:2:729,:) + 1i*h1(2:2:730,:); hh.off=h2(1:2:729,:) + 1i*h2(2:2:730,:);
vv.ref=v1(1:2:729,:) + 1i*v1(2:2:730,:); vv.off=v2(1:2:729,:) + 1i*v2(2:2:730,:);
xx.ref=x1(1:2:729,:) + 1i*x1(2:2:730,:); xx.off=x2(1:2:729,:) + 1i*x2(2:2:730,:);

hhref = hh.ref;hhoff = hh.off;
vvref = vv.ref;vvoff = vv.off;
xxref = xx.ref;xxoff = xx.off;

Pol_second_ground = [1;1;0];
Pol_second_ground = Pol_second_ground / norm(Pol_second_ground);
Pol_second_vegitation = [1;1;1];
Pol_second_vegitation = Pol_second_vegitation / norm(Pol_second_vegitation);

% Pol_fourth_ground = [1;1;0;-1;0;0;]; %ground
% Pol_fourth_vegitation = [1;1;1;-1;1;1;]; %vegitation

r = 7; c = 7;
L = r+c+1;
Lreshape = (r+c+1)^2;
[ylength,xlength] = size(hhoff);
g = zeros(ylength-r,xlength-c);
v = zeros(ylength-r,xlength-c);
n = zeros(ylength-r,xlength-c);
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = fliplr(r+1:ylength-r);
    for col = c+1:xlength-c;
        
        Var = hhref(row-r:row+r,col-c:col+c);
        s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = hhoff(row-r:row+r,col-c:col+c);
        s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = vvref(row-r:row+r,col-c:col+c);
        s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = vvoff(row-r:row+r,col-c:col+c);
        s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = xxref(row-r:row+r,col-c:col+c);
        s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = xxoff(row-r:row+r,col-c:col+c);
        s2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        
        R1 = s1*s1';  R2 = s1*s2';
        
        [eigenvec,eigenval] = eig(pinv(R1)*R2);
        
%         eigenveca = abs(eigenvec)
%         eigenvalaa = abs(angle(diag(eigenval)))

        [~,srt]=sort(abs(angle(diag(eigenval))),'descend');
        
        Leig_copol = abs(eigenvec(1,srt(1)))^2 + abs(eigenvec(2,srt(1)))^2;
        SLeig_copol = abs(eigenvec(1,srt(2)))^2 + abs(eigenvec(2,srt(2)))^2;
        
        n(row,col) = eigenval(srt(3),srt(3));
        
        if (Leig_copol >= SLeig_copol)
            
           g(row,col) = eigenval(srt(1),srt(1));
           v(row,col) = eigenval(srt(2),srt(2));
           
        else    
            
           g(row,col) = eigenval(srt(2),srt(2));
           v(row,col) = eigenval(srt(1),srt(1));
           
        end
        
    end
end

%% Plotting Results
figure(1); imagesc(angle(g)); title('angle(g)');
figure(2); imagesc(angle(v)); title('angle(v)');
figure(3); imagesc(angle(n)); title('angle(n)');
figure(5); imagesc(abs(g)); title('abs(g)');
figure(6); imagesc(abs(v)); title('abs(v)');
figure(7); imagesc(abs(n)); title('abs(n)');