%% Second order Algorithm Conditional
%% Initializations
clear;clc;

[h_data, v_data, x_data] = openSARdata();

row_window = 3; col_window = 3;

% L = row_window+col_window+1;

[ylength,xlength] = size(h_data.off);


phase.one = zeros(ylength,xlength);
phase.two = zeros(ylength,xlength);
phase.three = zeros(ylength,xlength);
phase.four = zeros(ylength,xlength);
phase.five = zeros(ylength,xlength);
phase.six = zeros(ylength,xlength);
phase.seven = zeros(ylength,xlength);
phase.eight = zeros(ylength,xlength);

Pol.one = [0000.1;0000.1;0000.1]/sqrt(0.0001);
Pol.two = [0;0;1]/sqrt(1);
Pol.three = [0;1;0]/sqrt(1);
Pol.four = [0;1;1]/sqrt(2);
Pol.five = [1;0;0]/sqrt(1);
Pol.six = [1;0;1]/sqrt(2);
Pol.seven = [1;1;0]/sqrt(2);
Pol.eight = [1;1;1]/sqrt(3);

% mask = ones(ylength,xlength)/(sqrt(ylength.^2 + xlength.^2));
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for column = 1:xlength;
        
        clear('R1','R2','S1','S2','A','uv','kk')
        
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
        
        polchg.one = abs(Pol.one'*u);
        polchg.two = abs(Pol.two'*u);
        polchg.three = abs(Pol.three'*u);
        polchg.four = abs(Pol.four'*u);
        polchg.five = abs(Pol.five'*u);
        polchg.six = abs(Pol.six'*u);
        polchg.seven = abs(Pol.seven'*u);
        polchg.eight = abs(Pol.eight'*u);
        
        [~,one] = sort(polchg.one,'descend');
        [~,two] = sort(polchg.two,'descend');
        [~,three] = sort(polchg.three,'descend');
        [~,four] = sort(polchg.four,'descend');
        [~,five] = sort(polchg.five,'descend');
        [~,six] = sort(polchg.six,'descend');
        [~,seven] = sort(polchg.seven,'descend');
        [~,eight] = sort(polchg.eight,'descend');
        
        
        phase.one(row,column) = (uv(one(1),one(1)));
        phase.two(row,column) = (uv(two(1),two(1)));
        phase.three(row,column) = (uv(three(1),three(1)));
        phase.four(row,column) = (uv(four(1),four(1)));
        phase.five(row,column) = (uv(five(1),five(1)));
        phase.six(row,column) = (uv(six(1),six(1)));
        phase.seven(row,column) = (uv(seven(1),seven(1)));
        phase.eight(row,column) = (uv(eight(1),eight(1)));
        
    end
end
%% Unwrapping Image
[unwrappedPhase.seven,unwrappedMag.seven] = QualityGuidedUnwrap2D( phase.seven);
[unwrappedPhase.eight,unwrappedMag.eight] = QualityGuidedUnwrap2D( phase.eight);
%% Plotting Results
HSIZE = 5;
SIGMA = .5;
% lowPass_phase.one = imfilter(angle(phase.one),fspecial('gaussian', HSIZE,SIGMA));
% lowPass_phase.two = imfilter(angle(phase.two),fspecial('gaussian', HSIZE,SIGMA));
% lowPass_phase.three = imfilter(angle(phase.three),fspecial('gaussian',HSIZE,SIGMA));
% lowPass_phase.four = imfilter(angle(phase.four),fspecial('gaussian', HSIZE,SIGMA));
% lowPass_phase.five = imfilter(angle(phase.five),fspecial('gaussian', HSIZE,SIGMA));
% lowPass_phase.six = imfilter(angle(phase.six),fspecial('gaussian', HSIZE,SIGMA));
lowPass_phase.seven = imfilter(unwrappedPhase.seven,fspecial('gaussian', HSIZE,SIGMA));
lowPass_phase.eight = imfilter(unwrappedPhase.eight,fspecial('gaussian', HSIZE,SIGMA));

% highPass_phase.one = angle(phase.one)-imfilter(angle(phase.one),fspecial('gaussian', HSIZE,SIGMA));
% highPass_phase.two = angle(phase.two)-imfilter(angle(phase.two),fspecial('gaussian', HSIZE,SIGMA));
% highPass_phase.three = angle(phase.three)-imfilter(angle(phase.three),fspecial('gaussian',HSIZE,SIGMA));
% highPass_phase.four = angle(phase.four)-imfilter(angle(phase.four),fspecial('gaussian', HSIZE,SIGMA));
% highPass_phase.five = angle(phase.five)-imfilter(angle(phase.five),fspecial('gaussian', HSIZE,SIGMA));
% highPass_phase.six = angle(phase.six)-imfilter(angle(phase.six),fspecial('gaussian', HSIZE,SIGMA));
highPass_phase.seven = angle(phase.seven)-imfilter(unwrappedPhase.seven,fspecial('gaussian', HSIZE,SIGMA));
highPass_phase.eight = angle(phase.eight)-imfilter(unwrappedPhase.eight,fspecial('gaussian', HSIZE,SIGMA));

% imtool(lowPass_phase.one);colormap(gray); %title('Pol.one = [0;0;0]/sqrt(0.00001);');
% imtool(lowPass_phase.two);colormap(gray); %title('Pol.two = [0;0;1]/sqrt(1);');
% imtool(lowPass_phase.three);colormap(gray); %title('Pol.three = [0;1;0]/sqrt(1);');
% imtool(lowPass_phase.four);colormap(gray); %title('Pol.four = [0;1;1]/sqrt(2);');
% imtool(lowPass_phase.five);colormap(gray); %title('Pol.five = [1;0;0]/sqrt(1);');
% imtool(lowPass_phase.six);colormap(gray); %title('Pol.six = [1;0;1]/sqrt(2);');
imtool(lowPass_phase.seven);colormap(gray);% title('Pol.seven = [1;1;0]/sqrt(2);');
imtool(lowPass_phase.eight);colormap(gray); %title('Pol.eight = [1;1;1]/sqrt(3);');


% imtool(highPass_phase.one);colormap(gray); %title('Pol.one = [0;0;0]/sqrt(0.00001);');
% imtool(highPass_phase.two);colormap(gray); %title('Pol.two = [0;0;1]/sqrt(1);');
% imtool(highPass_phase.three);colormap(gray); %title('Pol.three = [0;1;0]/sqrt(1);');
% imtool(highPass_phase.four);colormap(gray); %title('Pol.four = [0;1;1]/sqrt(2);');
% imtool(highPass_phase.five);colormap(gray); %title('Pol.five = [1;0;0]/sqrt(1);');
% imtool(highPass_phase.six);colormap(gray); %title('Pol.six = [1;0;1]/sqrt(2);');
imtool(highPass_phase.seven);colormap(gray);% title('Pol.seven = [1;1;0]/sqrt(2);');
imtool(highPass_phase.eight);colormap(gray); %title('Pol.eight = [1;1;1]/sqrt(3);');

