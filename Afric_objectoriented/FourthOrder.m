classdef FourthOrder < handle
    %FOURTHORDER Summary of this class
    %% Cumulant Algorithm conditional Fourth order
    
    %   Detailed explanation
    
    properties
        userinput_window_rowmax;
        userinput_window_columnmax;
  %  end
    %% Private properties of the ForthOrder Class
  %  properties (Access = private)
        ground_image;
        vegitation_image;
        noise_image;
        
        SMSAR_totalrows;
        SMSAR_totalcolums;
        
        SM_polinsar = struct();
    end
    %% Creating Contructor Function
    methods
        function obj = FourthOrder(ui_window_rowmax , ui_window_columnmax)
            if(ui_window_rowmax <= 0)
                error('Window Row must be greater than zereo')
            elseif (ui_window_columnmax <= 0)
                error('Window Column must be greater than zereo')
            end
            
            obj.userinput_window_rowmax = ui_window_rowmax;
            obj.userinput_window_columnmax = ui_window_columnmax;
        end
        %% Static Initial and plot calls Of the FourthOrder Class
        
        function retcode = plot_absolutevalue(obj)
            %% Plotting private image absolute values
            figure(1); imagesc(abs(obj.ground_image)); title('4D abs(g)');
            figure(2); imagesc(abs(obj.vegitation_image)); title('4D abs(v)');
            figure(3); imagesc(abs(obj.noise_image)); title('4D abs(n)');
            
            retcode = 0;
        end
        function retcode = plot_angle(obj)
            %% Plotting private image angles
            figure(4); imagesc(angle(obj.ground_image)); title('4D angle(g)');
            figure(5); imagesc(angle(obj.vegitation_image)); title('4D angle(v)');
            figure(6); imagesc(angle(obj.noise_image)); title('4D angle(n)');
            figure(7); imagesc(abs(obj.vegitation_image - abs(obj.ground_image))); title('4D abs(v)-abs(g)');
            figure(8); imagesc(angle(obj.vegitation_image) - angle(obj.ground_image)); title('4D angle(v)-angle(g)');
            
            retcode = 0;
        end
        function retcode = init(obj)
            %% initializing
            [obj.SM_polinsar.hh,obj.SM_polinsar.vv,obj.SM_polinsar.xx] =  openSARdata();      
            
            [obj.SMSAR_totalrows,obj.SMSAR_totalcolums] = size(obj.SM_polinsar.hh.off);
            
            obj.ground_image = zeros(obj.SMSAR_totalrows , obj.SMSAR_totalcolums);
            obj.vegitation_image = zeros(obj.SMSAR_totalrows, obj.SMSAR_totalcolums);
            obj.noise_image = zeros(obj.SMSAR_totalrows , obj.SMSAR_totalcolums);
            
            retcode = 0;
        end
        
        function retcode = run(obj)
            %% Espirit Martix Calculations Main Thread
            Row_info.window_rowmax = obj.userinput_window_rowmax;
            Col_info.window_colmax = obj.userinput_window_columnmax;
            
            for current_row = 1:obj.SMSAR_totalrows;
                for current_col = 1:obj.SMSAR_totalcolums;
                    
                    clear('R1','R2','S1','S2','uv','kk','s1','s2')
                    
                    Row_info.current_row = current_row;
                    Col_info.current_column = current_col;
                    
                    if(Row_info.current_row <= (Row_info.window_rowmax + 1))...
                            &&(Col_info.current_column <= ( Col_info.window_colmax + 1))    % condition 1
                        [s1,s2] = Average_Condition1(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row <= (Row_info.window_rowmax + 1))...
                            &&(Col_info.current_column > ( Col_info.window_colmax + 1))&&(Col_info.current_column < (obj.SMSAR_totalrows - ( Col_info.window_colmax - 1)))    % condition 2
                        [s1,s2] = Average_Condition2(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row <= (Row_info.window_rowmax + 1))...
                            &&(Col_info.current_column >= ((obj.SMSAR_totalcolums - ( Col_info.window_colmax - 1))))   % condition 3
                        [s1,s2] = Average_Condition3(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row  > (Row_info.window_rowmax + 1))&&(Row_info.current_row < (obj.SMSAR_totalrows - (Row_info.window_rowmax - 1)))...
                            &&(Col_info.current_column <= ( Col_info.window_colmax + 1))   % condition 4
                        [s1,s2] = Average_Condition4(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row > (Row_info.window_rowmax + 1))&&(Row_info.current_row < (obj.SMSAR_totalrows - (Row_info.window_rowmax -1)))...
                            &&(Col_info.current_column > ( Col_info.window_colmax + 1))&&(Col_info.current_column < (obj.SMSAR_totalcolums - ( Col_info.window_colmax-1)))  % condition 5
                        [s1,s2] = Average_Condition5(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row > (Row_info.window_rowmax + 1))&&(Row_info.current_row < (obj.SMSAR_totalrows - (Row_info.window_rowmax - 1)))...
                            &&(Col_info.current_column >= ((obj.SMSAR_totalcolums - ( Col_info.window_colmax - 1))))     % condition 6
                        [s1,s2] = Average_Condition6(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row >= (obj.SMSAR_totalrows - (Row_info.window_rowmax - 1)))...
                            &&(Col_info.current_column <= ( Col_info.window_colmax + 1))    % condition 7
                        [s1,s2] = Average_Condition7(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row >= (obj.SMSAR_totalrows - Row_info.window_rowmax))...
                            &&(Col_info.current_column > ( Col_info.window_colmax + 1))&&(Col_info.current_column < ((obj.SMSAR_totalcolums - ( Col_info.window_colmax - 1))))   % condition 8
                        [s1,s2] = Average_Condition8(Row_info,Col_info,obj.SM_polinsar);
                        
                    elseif(Row_info.current_row >= (obj.SMSAR_totalrows - (Row_info.window_rowmax - 1)))...
                            &&(Col_info.current_column >= ((obj.SMSAR_totalcolums - ( Col_info.window_colmax - 1))))   % condition 9
                        [s1,s2] = Average_Condition9(Row_info,Col_info,obj.SM_polinsar);
                    end
                    
                    S1 = [s1.h.*s1.h;
                        s1.h.*s1.v;
                        s1.h.*s1.x];
                    
                    S2 = [s2.h.*s2.h;
                        s2.h.*s2.v;
                        s2.h.*s2.x];
                    
                    R1 = S1*S1';
                    R2 = S1*S2';
                    
                    [~,uv] = eig(pinv(R1)*R2);
                    
                    [~,kk]=sort(angle(diag(uv)),'ascend');
                    
                    obj.ground_image(Row_info.current_row,Col_info.current_column) = uv(kk(1),kk(1));
                    obj.vegitation_image(Row_info.current_row,Col_info.current_column) = uv(kk(2),kk(2));
                    obj.noise_image(Row_info.current_row,Col_info.current_column) = uv(kk(3),kk(3));
                    
                end
            end
            retcode = 0;
            
            
        end
    end
end