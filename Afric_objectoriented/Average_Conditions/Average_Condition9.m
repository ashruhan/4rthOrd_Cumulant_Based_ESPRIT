%% Condition 9
function [S1,S2] = Average_Condition9(Row_info,Col_info,SM_polinsar)

SM_polinsar.hh.ref = SM_polinsar.hh.ref((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);

[row_total,col_total] = size(SM_polinsar.hh.ref); 
window_area = row_total*col_total;

S1.h = reshape(SM_polinsar.hh.ref,1,window_area);

SM_polinsar.hh.off = SM_polinsar.hh.off((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);
S2.h = reshape(SM_polinsar.hh.off,1,window_area);

SM_polinsar.vv.ref = SM_polinsar.vv.ref((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);
S1.v = reshape(SM_polinsar.vv.ref,1,window_area);

SM_polinsar.vv.off = SM_polinsar.vv.off((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);
S2.v = reshape(SM_polinsar.vv.off,1,window_area);

SM_polinsar.xx.ref = SM_polinsar.xx.ref((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);
S1.x = reshape(SM_polinsar.xx.ref,1,window_area);

SM_polinsar.xx.off = SM_polinsar.xx.off((Row_info.current_row - Row_info.window_rowmax):end,...
    (Col_info.current_column - Col_info.window_colmax):end);
S2.x = reshape(SM_polinsar.xx.off,1,window_area);
end
