%% Cumulant Algorithm conditional Fourth order

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
        
        S1 = S1*exp(0.5);
        S2 = S2*exp(0.5);
        
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
figure(1); imagesc(angle(ground)); title('4D angle(g)');
figure(2); imagesc(angle(vegitation)); title('4D angle(v)');
figure(3); imagesc(angle(noise)); title('4D angle(n)');
figure(4); imagesc(abs(ground)); title('4D abs(g)');
figure(5); imagesc(abs(vegitation)); title('4D abs(v)');
figure(6); imagesc(abs(noise)); title('4D abs(n)');
figure(7); imagesc(abs(vegitation - abs(ground))); title('4D abs(v)-abs(g))');
figure(8); imagesc(angle(vegitation) - angle(ground)); title('4D angle(v)-angle(g)');