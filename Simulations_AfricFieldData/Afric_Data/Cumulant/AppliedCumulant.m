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

X = zeros(1,xlength);
Y = zeros(1,ylength);

halfrow = 7; halfcol = 7;
base = halfrow+halfcol+1;
Window_optimal = (base)^2;
eye_4 = -0.15;

ground = zeros(ylength-halfrow,xlength-halfcol);
vegetation = zeros(ylength-halfrow,xlength-halfcol);
%% Cumulant Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = halfrow+1:ylength-halfrow;
    for col = halfcol+1:xlength-halfcol;
        
        Var(1:base,1:base) = hhref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sh1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = hhoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sh2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = vvref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sv1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = vvoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sv2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = xxref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sx1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = xxoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sx2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        s1_Noise(1,:) = sh1;
        s1_Noise(2,:) = sv1;
        s1_Noise(3,:) = sx1;
        
        s2_Noise(1,:) = sh2;
        s2_Noise(2,:) = sv2;
        s2_Noise(3,:) = sx2;
        
        [ Cumulant_11,Cumulant_12] = Cumulant( s1_Noise,s2_Noise ,Window_optimal);      
        
        [eigenvec_4,eigenval_4] = eig(pinv(Cumulant_11 + eye_4*eye(6))*Cumulant_12);
        
        [~,srt_4] = sort(abs(diag(eigenval_4)),'descend');
        
        valord = 1;
        vecord = 1;
        
        LeigTemp  = (abs(eigenval_4(srt_4(1),srt_4(1))))^valord....
            *(abs(eigenvec_4(3,srt_4(1)))^vecord....
            + abs(eigenvec_4(5,srt_4(1)))^vecord....
            + abs(eigenvec_4(6,srt_4(1)))^vecord);
        
        SLeigTemp = (abs(eigenval_4(srt_4(2),srt_4(2))))^valord....
            *(abs(eigenvec_4(3,srt_4(2)))^vecord....
            + abs(eigenvec_4(5,srt_4(2)))^vecord....
            + abs(eigenvec_4(6,srt_4(2)))^vecord);
        
        if LeigTemp >= SLeigTemp
            
            vegetation(row,col) = eigenval_4(srt_4(1),srt_4(1));
            
            ground(row,col) = eigenval_4(srt_4(2),srt_4(2));
             
        else
            
            vegetation(row,col) = eigenval_4(srt_4(2),srt_4(2));
            ground(row,col) = eigenval_4(srt_4(1),srt_4(1));
        end       
    end
end
%% Plotting gv Results
gndscl = size(unique(reshape(ground,size(ground,1)*size(ground,2),size(ground,3))),1);
vegscl = size(unique(reshape(vegetation,size(ground,1)*size(ground,2),size(ground,3))),1);

figure(1); imshow(0.5*angle(ground),'DisplayRange',[-pi pi],'Colormap',jet(gndscl)); title('4th Ord angle(g)');
figure(2); imshow(0.5*angle(vegetation),'DisplayRange',[-pi pi],'Colormap',jet(vegscl)); title('4th Ord angle(v)');
figure(3); imshow(sqrt(abs(ground)),'DisplayRange',[0 mean(mean(sqrt(abs(ground))))],'Colormap',jet(gndscl)); title('4th Ord abs(g)');
figure(4); imshow(sqrt(abs(vegetation)),'DisplayRange',[0 mean(mean(sqrt(abs(vegetation))))],'Colormap',jet(vegscl)); title('4th Ord abs(v)');

% figure(1); imagesc(0.5*angle(ground),[0 2*pi]); title('4th Ord angle(g)');
% figure(2); imagesc(0.5*angle(vegetation),[0 2*pi]); title('4th Ord angle(v)');
% figure(3); imagesc(sqrt(abs(ground)),[0 10]); title('4th Ord abs(g)');
% figure(4); imagesc(sqrt(abs(vegetation)),[0 10]); title('4th Ord abs(v)');