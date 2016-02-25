%% Espirit Algorithm conditional Second order

%    window_col       window_col
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #   window_row
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #
%   $ $ $ $ $ $ $ * $ $ $ $ $ $ $
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #                                                      
%   # # # # # # # $ # # # # # # #   window_row
%   # # # # # # # $ # # # # # # #                                                                 
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #
%   # # # # # # # $ # # # # # # #



%% Initializations
[hh,vv,xx] =  openSARdata();
window_row = 7; 
window_col = 7;

[total_rows,total_cols] = size(hh.off);

ground = zeros(total_rows,total_cols);
vegitation = zeros(total_rows,total_cols);
noise = zeros(total_rows,total_cols);
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for current_row = 1:total_rows;
    for current_col = 1:total_cols;
        
        clear('R1','R2','S1','S2','uv','kk')
        
        Row_info.window_row = window_row;
        Row_info.row = current_row;
        
        Col_info.window_col = window_col;
        Col_info.col = current_col;
        
        if(current_row <= (window_row + 1))...
                &&(current_col <= (window_col + 1))     % condition 1
            [s1,s2] = Average_Condition1(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row <= (window_row + 1))...
                &&(current_col > (window_col + 1))&&(current_col < (total_cols - (window_col - 1)))    % condition 2
            [s1,s2] = Average_Condition2(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row <= (window_row + 1))...
                &&(current_col >= (total_cols - (window_col - 1)))   % condition 3
            [s1,s2] = Average_Condition3(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row  > (window_row + 1))&&(current_row < (total_rows - (window_row - 1)))...
                &&(current_col <= (window_col + 1))   % condition 4
            [s1,s2] = Average_Condition4(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row > (window_row + 1))&&(current_row < total_rows - (window_row -1))...
            &&(current_col > (window_col + 1))&&(current_col < total_cols - (window_col-1))  % condition 5
            [s1,s2] = Average_Condition5(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row > (window_row + 1))&&(current_row < (total_rows - window_row))...
            &&(current_col >= (total_cols - (window_col - 1)))     % condition 6
            [s1,s2] = Average_Condition6(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row >= (total_rows - (window_row - 1)))...
                &&(current_col <= (window_col + 1))    % condition 7
            [s1,s2] = Average_Condition7(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row >= (total_rows - window_row))...
                &&(current_col > (window_col + 1))&&(current_col < (total_cols - (window_col - 1)))   % condition 8
            [s1,s2] = Average_Condition8(Row_info,Col_info,hh,vv,xx);
            
        elseif(current_row >= (total_rows - (window_row - 1)))...
                &&(current_col >= (total_cols - (window_col - 1)))   % condition 9
            [s1,s2] = Average_Condition9(Row_info,Col_info,hh,vv,xx);
        end
        
        S1(1,:) = s1.h; S2(1,:) = s2.h;
        S1(2,:) = s1.v; S2(2,:) = s2.v;
        S1(3,:) = s1.x; S2(3,:) = s2.x;
        
        
        R1 = S1*S1';  
        
        R2 = S1*S2';
        
        [~,uv] = eig(pinv(R1)*R2);
        
        [~,kk]=sort(angle(diag(uv)),'ascend');
        
        ground(current_row,current_col) = uv(kk(1),kk(1));
        vegitation(current_row,current_col) = uv(kk(2),kk(2));
        noise(current_row,current_col) = uv(kk(3),kk(3));
        
    end
end

%% Plotting Results
 
figure(9); imagesc(angle(ground)); title('2D angle(g)');
figure(10); imagesc(angle(vegitation)); title('2D angle(v)');
figure(11); imagesc(angle(noise)); title('2D angle(n)');
figure(12); imagesc(abs(ground)); title('2D abs(g)');
figure(13); imagesc(abs(vegitation)); title('2D abs(v)');
figure(14); imagesc(abs(noise)); title('2D abs(n)');
figure(15); imagesc(abs(vegitation - abs(ground))); title('2D abs(v)-abs(g)');
figure(16); imagesc(angle(vegitation) - angle(ground)); title('2D angle(v)-angle(g)');