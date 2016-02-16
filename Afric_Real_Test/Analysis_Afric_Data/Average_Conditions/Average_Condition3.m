%% Condition 3
function [S1,S2] = Average_Condition3(Row_info,Col_info,hh,vv,xx)

hh.ref_temp = hh.ref(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);

[rows_total,col_total] = size(hh.ref_temp); 

window_area = rows_total*col_total;
S1.h = reshape(hh.ref_temp,1,window_area);

hh.ref_temp = hh.off(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);
S2.h = reshape(hh.ref_temp,1,window_area);

vv.off_temp = vv.ref(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);
S1.v = reshape(vv.off_temp,1,window_area);

vv.off_temp = vv.off(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);
S2.v = reshape(vv.off_temp ,1,window_area);

xx.ref_temp = xx.ref(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);
S1.x = reshape(xx.ref_temp,1,window_area);

xx.off_temp = xx.off(1:(Row_info.row + Row_info.window_row),(Col_info.col - Col_info.window_col):end);
S2.x = reshape(xx.off_temp,1,window_area);

end