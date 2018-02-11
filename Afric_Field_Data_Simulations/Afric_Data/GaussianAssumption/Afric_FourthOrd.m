%%Espirit Algorithm Fourth order
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

[ylength,xlength] = size(hhoff);
r = 7; c = 7;
L = r+c+1;
Lreshape = (r+c+1)^2;

ground_4 = zeros(ylength-r,xlength-c);
vegetation_4 = zeros(ylength-r,xlength-c);
%% Cumulant Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = r+1:ylength-r;
    for col = c+1:xlength-c;
        
        Var(1:L,1:L) = hhref(row-r:row+r,col-c:col+c);
        sh1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = hhoff(row-r:row+r,col-c:col+c);
        sh2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = vvref(row-r:row+r,col-c:col+c);
        sv1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = vvoff(row-r:row+r,col-c:col+c);
        sv2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = xxref(row-r:row+r,col-c:col+c);
        sx1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = xxoff(row-r:row+r,col-c:col+c);
        sx2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        S1 = [sh1.*sh1;sv1.*sv1;sx1.*sx1;sh1.*sv1;sh1.*sx1;sv1.*sx1;];
        S2 = [sh2.*sh2;sv2.*sv2;sx2.*sx2;sh2.*sv2;sh2.*sx2;sv2.*sx2;];
        
        R1_4 = S1*S1';
        
        R2_4 = S1*S2';
        
        [~,eigenvalCov11_4] = eig(R1_4,'nobalance');
        
        [~,srtCov_4] = sort(abs(diag(eigenvalCov11_4)),'descend');
        
        eye_4 = mean((diag(eigenvalCov11_4)))/max(abs(diag(eigenvalCov11_4)))^2;
        
        [eigenvec_4,eigenval_4] = eig((pinv(R1_4 + eye_4*eye(6)))*R2_4,'nobalance');
        
        [~,srt_abs_4] = sort(abs(diag(eigenval_4)),'descend');
        
        [~,srt_angle_4] = sort(abs(angle(diag(eigenval_4(srt_abs_4(1:3),srt_abs_4(1:3))))),'descend');
        
        LeigTemp  = (abs(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))^2....
            *(abs(eigenvec_4(3,srt_abs_4(srt_angle_4(1))))^2....
            + abs(eigenvec_4(5,(srt_angle_4(1))))^2....
            + abs(eigenvec_4(6,srt_abs_4(srt_angle_4(1))))^2);
        
        SLeigTemp = (abs(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))^2....
            *(abs(eigenvec_4(3,srt_abs_4(srt_angle_4(3))))^2....
            + abs(eigenvec_4(5,srt_abs_4(srt_angle_4(3))))^2....
            + abs(eigenvec_4(6,srt_abs_4(srt_angle_4(3))))^2);
        
        if ((LeigTemp >= SLeigTemp)&&(abs(eigenvec_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))...
                >= abs(eigenvec_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))) )
            
            ground_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)));
            
            vegetation_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)));
            
        elseif((LeigTemp >= SLeigTemp)&&(abs(eigenvec_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))...
                >= abs(eigenvec_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))) )
            
            ground_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)));
            
            vegetation_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)));
            
            
        elseif (SLeigTemp >= LeigTemp)
            
            ground_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)));
            
            vegetation_4(row,col) = eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)));
            
        end
        
    end
end
%% Plotting gv Results
figure(1); imagesc(angle(ground_4)); title('4th Order Ground Interferometric Phase');colorbar;
figure(2); imagesc(angle(vegetation_4)); title('4th Order Vegetation Interferometric Phase');colorbar;
figure(3); imagesc(angle(vegetation_4) - angle(ground_4)); title('4th Order V-G');colorbar;
figure(4); imagesc(abs(ground_4)); title('4th Order Ground Magnitude');colorbar;
figure(5); imagesc(abs(vegetation_4)); title('4th Order Vegetation Magnitude');colorbar;