%% Espirit Algorithm conditional Second and Fourth Order
%% Initializations
clc;clear;

[hh,vv,xx] =  openSARdata();
row_window = 10; col_window = 10;
[ylength,xlength] = size(hh.off);

inter.h.complex = zeros(ylength,xlength);
inter.v.complex = zeros(ylength,xlength);
inter.x.complex = zeros(ylength,xlength);

secondOrd.largest.complex = zeros(ylength,xlength);
secondOrd.medium.complex = zeros(ylength,xlength);

fourthOrd.largest.complex = zeros(ylength,xlength);
fourthOrd.medium.complex = zeros(ylength,xlength);

%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for column = 1:xlength;
        
        clear('R1','R2','S1','S2','A','uv','kk')
        R.r = row_window;R.row = row;R.ylength = ylength;
        C.c = col_window;C.col = column;C.xlength = xlength;
        
        if(row<row_window+1)&&(column<col_window+1)     % condition 1
            [s1,s2] = Average_Condition1(R,C,hh,vv,xx);
            
        elseif(row<row_window+1)&&(col_window+1<column)&&(column<xlength-col_window)    % condition 2
            [s1,s2] = Average_Condition2(R,C,hh,vv,xx);
            
        elseif(row<row_window+1)&&(column>xlength-col_window)   % condition 3
            [s1,s2] = Average_Condition3(R,C,hh,vv,xx);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(column<col_window+1)   % condition 4
            [s1,s2] = Average_Condition4(R,C,hh,vv,xx);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(col_window+1<column)&&(column<xlength-col_window)  % condition 5
            [s1,s2] = Average_Condition5(R,C,hh,vv,xx);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(column>xlength-col_window)     % condition 6
            [s1,s2] = Average_Condition6(R,C,hh,vv,xx);
            
        elseif(row>ylength-row_window)&&(column<col_window+1)    % condition 7
            [s1,s2] = Average_Condition7(R,C,hh,vv,xx);
            
        elseif(row>ylength-row_window)&&(col_window+1<column)&&(column<xlength-col_window)   % condition 8
            [s1,s2] = Average_Condition8(R,C,hh,vv,xx);
            
        elseif(row>ylength-row_window)&&(column>xlength-col_window)   % condition 9
            [s1,s2] = Average_Condition9(R,C,hh,vv,xx);
        end
        %% Interferometry
        inter.h.complex(row,column) = s1.h*s2.h'/sqrt((s1.h*s1.h')*(s2.h*s2.h'));
        inter.v.complex(row,column) = s1.v*s2.v'/sqrt((s1.v*s1.v')*(s2.v*s2.v'));
        inter.x.complex(row,column) = s1.x*s2.x'/sqrt((s1.x*s1.x')*(s2.x*s2.x'));
        
        %% Second Order ESPIRIT
        S1(1,:) = s1.h;
        S1(2,:) = s1.v;
        S1(3,:) = s1.x;
        
        S2(1,:) = s2.h;
        S2(2,:) = s2.v;
        S2(3,:) = s2.x;
        
        R1 = S1*S1';
        R2 = S1*S2';
        A = pinv(R1)*R2;
        
        [~,uv] = eig(A);
        
        [~,kk] = sort(abs(diag(uv)),'ascend');
        
        secondOrd.largest.complex(row,column) = uv(kk(1),kk(1));
        secondOrd.medium.complex(row,column) = uv(kk(2),kk(2));
        
        %% Fourth Order ESPIRIT
        clear('R1','R2','S1','S2','A','uv','kk')
        S1 = [s1.h.*s1.h;
            s1.v.*s1.v;
            s1.x.*s1.x;
            s1.h.*s1.v;
            s1.h.*s1.x;
            s1.v.*s1.x];
        
        S2 = [s2.h.*s2.h;
            s2.v.*s2.v;
            s2.x.*s2.x;
            s2.h.*s2.v;
            s2.h.*s2.x;
            s2.v.*s2.x];
        
        
        R1 = S1*S1';
        R2 = S1*S2';
        
        A = pinv(R1)*R2;
        
        [~,uv] = eig(A);
        
        [~,kk] = sort(abs(diag(uv)),'ascend');
        
        fourthOrd.largest.complex(row,column) = uv(kk(1),kk(1));
        fourthOrd.medium.complex(row,column) = uv(kk(2),kk(2));
        
    end
end
%% Returning Absolute and Phase Results
secondOrd.largest.phase.phase = angle(secondOrd.largest.complex);
secondOrd.medium.phase.phase = angle(secondOrd.medium.complex);
fourthOrd.largest.phase.phase = angle(fourthOrd.largest.complex);
fourthOrd.medium.phase.phase = angle(fourthOrd.medium.complex);

secondOrd.largest.abs.abs = abs(secondOrd.largest.complex);
secondOrd.medium.abs.abs = abs(secondOrd.medium.complex);
fourthOrd.largest.abs.abs = abs(fourthOrd.largest.complex);
fourthOrd.medium.abs.abs = abs(fourthOrd.medium.complex);

inter.h.phase.phase = angle(inter.h.complex);
inter.v.phase.phase = angle(inter.v.complex);
inter.x.phase.phase = angle(inter.x.complex);
inter.h.abs.abs = abs(inter.h);
inter.v.abs.abs = abs(inter.v);
inter.x.abs.abs = abs(inter.x);
%% Unwrapping Image
% [unwrappedPhase.seven,unwrappedMag.seven] = QualityGuidedUnwrap2D( phase.seven);
% [unwrappedPhase.eight,unwrappedMag.eight] = QualityGuidedUnwrap2D( phase.eight);
%% Low Pass and High Pass Filtering Results
HSIZE = 10;
SIGMA = .5;

inter.h.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.h.abs.abs);
inter.v.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.v.abs.abs);
inter.x.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.x.abs.abs);

inter.h.abs.highpass = inter.h.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.h.abs.abs);
inter.v.abs.highpass = inter.v.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.v.abs.abs);
inter.x.abs.highpass = inter.x.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.x.abs.abs);

secondOrd.largest.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.abs.abs);
secondOrd.medium.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.abs.abs);
fourthOrd.largest.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.abs.abs);
fourthOrd.medium.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.abs.abs);

secondOrd.largest.abs.highpass = secondOrd.largest.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.abs.abs);
secondOrd.medium.abs.highpass = secondOrd.medium.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.abs.abs);
fourthOrd.largest.abs.highpass = fourthOrd.largest.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.abs.abs);
fourthOrd.medium.abs.highpass =  fourthOrd.medium.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.abs.abs);


inter.h.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.h.phase.phase);
inter.v.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.v.phase.phase);
inter.x.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),inter.x.phase.phase);

inter.h.phase.highpass = inter.h.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.h.phase.phase);
inter.v.phase.highpass = inter.v.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.v.phase.phase);
inter.x.phase.highpass = inter.x.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),inter.x.phase.phase);

secondOrd.largest.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.phase.phase);
secondOrd.medium.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.phase.phase);
fourthOrd.largest.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.phase.phase);
fourthOrd.medium.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.phase.phase);

secondOrd.largest.phase.highpass = secondOrd.largest.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.phase.phase);
secondOrd.medium.phase.highpass = secondOrd.medium.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.phase.phase);
fourthOrd.largest.phase.highpass = fourthOrd.largest.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.phase.phase);
fourthOrd.medium.phase.highpass = fourthOrd.medium.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.phase.phase);

%% Finding Errors in ESPIRIT Using Interferometry Results
error.second.largest.h = secondOrd.largest.phase.lowpass - inter.h.phase.lowpass;
error.second.largest.v = secondOrd.largest.phase.lowpass - inter.v.phase.lowpass;
error.second.largest.x = secondOrd.largest.phase.lowpass - inter.x.phase.lowpass;

error.second.medium.h = secondOrd.medium.phase.lowpass - inter.h.phase.lowpass;
error.second.medium.v = secondOrd.medium.phase.lowpass - inter.v.phase.lowpass;
error.second.medium.x = secondOrd.medium.phase.lowpass - inter.x.phase.lowpass;


%% Plotting Results

figure(1);imshow(error.second.largest.h);colormap(parula);title('error.second.largest.h');
figure(1);imshow(error.second.largest.v);colormap(parula);title('error.second.largest.v');
figure(1);imshow(error.second.largest.x);colormap(parula);title('error.second.largest.x');
figure(1);imshow(error.second.medium.h);colormap(parula);title('error.second.medium.h');
figure(1);imshow(error.second.medium.v);colormap(parula);title('error.second.medium.v');
figure(1);imshow(error.second.medium.x);colormap(parula);title('error.second.medium.x');