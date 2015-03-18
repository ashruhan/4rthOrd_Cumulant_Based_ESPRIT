%% Second order Algorithm Conditional
%% Initializations
clear;clc;

[h_data, v_data, x_data] = openSARdata();

row_window = 10; col_window = 10;

[ylength,xlength] = size(h_data.off);

inter.h = zeros(ylength,xlength);
inter.v = zeros(ylength,xlength);
inter.x = zeros(ylength,xlength);

%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for column = 1:xlength;
        
%         clear('s1','s2')
        
        R.r = row_window; R.row = row; R.ylength = ylength;
        C.c = col_window; C.col = column; C.xlength = xlength;
        
        if(row<row_window+1)&&(column<col_window+1)     % condition 1
            [s1,s2] = Average_Condition1(R,C,h_data,v_data,x_data);
            
        elseif(row<row_window+1)&&(col_window+1<column)&&(column<xlength-col_window)    % condition 2
            [s1,s2] = Average_Condition2(R,C,h_data,v_data,x_data);
            
        elseif(row<row_window+1)&&(column>xlength-col_window)   % condition 3
            [s1,s2] = Average_Condition3(R,C,h_data,v_data,x_data);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(column<col_window+1)   % condition 4
            [s1,s2] = Average_Condition4(R,C,h_data,v_data,x_data);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(col_window+1<column)&&(column<xlength-col_window)  % condition 5
            [s1,s2] = Average_Condition5(R,C,h_data,v_data,x_data);
            
        elseif(row_window+1<row)&&(row<ylength-row_window)&&(column>xlength-col_window)     % condition 6
            [s1,s2] = Average_Condition6(R,C,h_data,v_data,x_data);
            
        elseif(row>ylength-row_window)&&(column<col_window+1)    % condition 7
            [s1,s2] = Average_Condition7(R,C,h_data,v_data,x_data);
            
        elseif(row>ylength-row_window)&&(col_window+1<column)&&(column<xlength-col_window)   % condition 8
            [s1,s2] = Average_Condition8(R,C,h_data,v_data,x_data);
            
        elseif(row>ylength-row_window)&&(column>xlength-col_window)   % condition 9
            [s1,s2] = Average_Condition9(R,C,h_data,v_data,x_data);
        end
        
        inter.h(row,column) = s1.h*s2.h'/sqrt((s1.h*s1.h')*(s2.h*s2.h'));
        inter.v(row,column) = s1.v*s2.v'/sqrt((s1.v*s1.v')*(s2.v*s2.v'));
        inter.x(row,column) = s1.x*s2.x'/sqrt((s1.x*s1.x')*(s2.x*s2.x'));
        
    end
end
phase_inter.h = angle(inter.h);
phase_inter.v = angle(inter.v);
phase_inter.x = angle(inter.x);
abs_inter.h = abs(inter.h);
abs_inter.v = abs(inter.v);
abs_inter.x = abs(inter.x);
%% Unwrapping Image
% [unwrappedPhase.seven,unwrappedMag.seven] = QualityGuidedUnwrap2D( phase.seven);
% [unwrappedPhase.eight,unwrappedMag.eight] = QualityGuidedUnwrap2D( phase.eight);
%% Plotting Results
HSIZE = 10;
SIGMA = .5;

lowPass_INTER_abs.h = filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.h);
lowPass_INTER_abs.v = filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.v);
lowPass_INTER_abs.x = filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.x);

highPass_INTER_abs.h = abs_inter.h-filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.h);
highPass_INTER_abs.v = abs_inter.v-filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.v);
highPass_INTER_abs.x = abs_inter.x-filter2(fspecial('gaussian', HSIZE,SIGMA),abs_inter.x);

lowPass_phase.h = filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.h);
lowPass_phase.v = filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.v);
lowPass_phase.x = filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.x);

highPass_phase.h = phase_inter.h-filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.h);
highPass_phase.v = phase_inter.v-filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.v);
highPass_phase.x = phase_inter.x-filter2(fspecial('gaussian', HSIZE,SIGMA),phase_inter.x);

figure(1);imshow(lowPass_INTER_abs.h);colormap(parula);title('lowPass.abs.h');
figure(2);imshow(lowPass_INTER_abs.v);colormap(parula);title('lowPass.abs.v');
figure(3);imshow(lowPass_INTER_abs.x);colormap(parula);title('lowPass.abs.x');

figure(4);imshow(highPass_INTER_abs.h);colormap(parula);title('highPass.abs.h');
figure(5);imshow(highPass_INTER_abs.v);colormap(parula);title('highPass.abs.v');
figure(6);imshow(highPass_INTER_abs.x);colormap(parula);title('highPass.abs.x');

figure(7);imshow(lowPass_phase.h);colormap(parula);title('lowPass.phase.h');
figure(8);imshow(lowPass_phase.v);colormap(parula);title('lowPass.phase.v');
figure(9);imshow(lowPass_phase.x);colormap(parula);title('lowPass.phase.x');

figure(10);imshow(highPass_phase.h);colormap(parula);title('highPass.phase.h');
figure(11);imshow(highPass_phase.v);colormap(parula);title('highPass.phase.v');
figure(12);imshow(highPass_phase.x);colormap(parula);title('highPass.phase.x');