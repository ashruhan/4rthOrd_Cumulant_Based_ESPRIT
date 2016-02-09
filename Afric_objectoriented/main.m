%% Main matlab script Executing High Level Classes
% This script will be used to call second or fourth order classes
% to (run) each class, a window column and window row needs be be
% initialized as seen below.

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

%% Initialize Window Row and Window Column
% If comparing Second and Fourth Use the same dimensions
clear; clc;
window_row = 7;
window_column = 7;

%% Initialize classes
FO = FourthOrder(window_row,window_column);

retcode = init(FO);
if (retcode)
    error('Error Initializing FourthOrder Init');
end

%% Execute main thread
retcode = run(FO);
if (retcode)
    error('Error Running Main thread of FourthOrder');
end

%% Imaging results
retcode = plot_absolutevalue(FO);
if (retcode)
    error('Error Imaging absolute values');
end

retcode = plot_angle(FO);
if (retcode)
    error('Error Imaging angles');
end

