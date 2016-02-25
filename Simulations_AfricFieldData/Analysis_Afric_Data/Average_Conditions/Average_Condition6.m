%% Condition 6
function [S1,S2] = Average_Condition6(Row_info,Col_info,hh,vv,xx)

hh.ref_temp = hh.ref((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);

[row_total,col_total] = size(hh.ref_temp); 

window_area = row_total*col_total;

S1.h = reshape(hh.ref_temp,1,window_area);

hh.off_temp = hh.off((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);
S2.h = reshape(hh.off_temp,1,window_area);

vv.ref_temp = vv.ref((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);
S1.v = reshape(vv.ref_temp,1,window_area);

vv.off_temp = vv.off((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);
S2.v = reshape(vv.off_temp,1,window_area);

xx.ref_temp = xx.ref((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);
S1.x = reshape(xx.ref_temp,1,window_area);

xx.off_temp = xx.off((Row_info.row - Row_info.window_row):(Row_info.row + Row_info.window_row)...
    ,(Col_info.col - Col_info.window_col):end);
S2.x = reshape(xx.off_temp,1,window_area);

end