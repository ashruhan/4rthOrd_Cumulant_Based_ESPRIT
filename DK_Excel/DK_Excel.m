%% Initializations
clear;clc;
FILE_NAME = 'Example.csv';

%% Make sure the file is there and name is correct

if exist(FILE_NAME,'file')
    temp_array = xlsread(FILE_NAME);
else
    error('Make sure your excel file is in the same directory as this .m file');
end

%% Extract meaningfull data
data.open_bell = temp_array(:,1);
data.high = temp_array(:,2);
data.low = temp_array(:,3);
data.close_bell = temp_array(:,4);
data.volume = temp_array(:,5);
data.adj_close = temp_array(:,6);

%% Market Trend Data

% Your code here;
