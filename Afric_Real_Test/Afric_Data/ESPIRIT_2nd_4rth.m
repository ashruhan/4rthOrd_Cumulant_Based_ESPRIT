%% Espirit Algorithm conditional Second and Fourth Order
%% Initializations
clc;clear;

[hh,vv,xx] =  openSARdata();
row_window = 10; col_window = 10;

[ylength,xlength] = size(hh.off);

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

secondOrd.largest.phase.phase = angle(secondOrd.largest.complex);
secondOrd.medium.phase.phase = angle(secondOrd.medium.complex);
fourthOrd.largest.phase.phase = angle(fourthOrd.largest.complex);
fourthOrd.medium.phase.phase = angle(fourthOrd.medium.complex);

secondOrd.largest.abs.abs = abs(secondOrd.largest.complex);
secondOrd.medium.abs.abs = abs(secondOrd.medium.complex);
fourthOrd.largest.abs.abs = abs(fourthOrd.largest.complex);
fourthOrd.medium.abs.abs = abs(fourthOrd.medium.complex);
%% Unwrapping Image
% [unwrappedPhase.seven,unwrappedMag.seven] = QualityGuidedUnwrap2D( phase.seven);
% [unwrappedPhase.eight,unwrappedMag.eight] = QualityGuidedUnwrap2D( phase.eight);
%% Plotting Results
HSIZE = 10;
SIGMA = .5;

secondOrd.largest.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.abs.abs);
secondOrd.medium.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.abs.abs);
fourthOrd.largest.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.abs.abs);
fourthOrd.medium.abs.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.abs.abs);

secondOrd.largest.abs.highpass = secondOrd.largest.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.abs.abs);
secondOrd.medium.abs.highpass = secondOrd.medium.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.abs.abs);
fourthOrd.largest.abs.highpass = fourthOrd.largest.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.abs.abs);
fourthOrd.medium.abs.highpass =  fourthOrd.medium.abs.abs - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.abs.abs);

secondOrd.largest.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.phase.phase);
secondOrd.medium.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.phase.phase);
fourthOrd.largest.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.phase.phase);
fourthOrd.medium.phase.lowpass = filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.phase.phase);

secondOrd.largest.phase.highpass = secondOrd.largest.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.largest.phase.phase);
secondOrd.medium.phase.highpass = secondOrd.medium.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),secondOrd.medium.phase.phase);
fourthOrd.largest.phase.highpass = fourthOrd.largest.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.largest.phase.phase);
fourthOrd.medium.phase.highpass = fourthOrd.medium.phase.phase - filter2(fspecial('gaussian', HSIZE,SIGMA),fourthOrd.medium.phase.phase);

figure(1);imshow(secondOrd.largest.abs.lowpass);colormap(parula);title('secondOrd.largest.abs.lowpass');
figure(2);imshow(secondOrd.medium.abs.lowpass);colormap(parula);title('secondOrd.medium.abs.lowpass');
figure(3);imshow(fourthOrd.largest.abs.lowpass);colormap(parula);title('fourthOrd.largest.abs.lowpass');
figure(4);imshow(fourthOrd.medium.abs.lowpass);colormap(parula);title('fourthOrd.medium.abs.lowpass');

figure(5);imshow(secondOrd.largest.abs.highpass);colormap(parula);title('secondOrd.largest.abs.highpass');
figure(6);imshow(secondOrd.medium.abs.highpass);colormap(parula);title('secondOrd.medium.abs.highpass');
figure(7);imshow(fourthOrd.largest.abs.highpass);colormap(parula);title('fourthOrd.largest.abs.highpass');
figure(8);imshow(fourthOrd.medium.abs.highpass);colormap(parula);title('fourthOrd.medium.abs.highpass');

figure(9);imshow(secondOrd.largest.phase.lowpass);colormap(parula);title('secondOrd.largest.phase.lowpass');
figure(10);imshow(secondOrd.medium.phase.lowpass);colormap(parula);title('secondOrd.medium.phase.lowpass');
figure(11);imshow(fourthOrd.largest.phase.lowpass);colormap(parula);title('fourthOrd.largest.phase.lowpass');
figure(12);imshow(fourthOrd.medium.phase.lowpass);colormap(parula);title('fourthOrd.medium.phase.lowpass');

figure(13);imshow(secondOrd.largest.phase.highpass);colormap(parula);title('secondOrd.largest.phase.highpass');
figure(14);imshow(secondOrd.medium.phase.highpass);colormap(parula);title('secondOrd.medium.phase.highpass');
figure(15);imshow(fourthOrd.largest.phase.highpass);colormap(parula);title('fourthOrd.largest.phase.highpass');
figure(16);imshow(fourthOrd.medium.phase.highpass);colormap(parula);title('fourthOrd.medium.phase.highpass');

