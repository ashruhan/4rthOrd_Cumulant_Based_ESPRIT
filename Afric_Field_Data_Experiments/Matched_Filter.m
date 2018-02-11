function Matched_Filter

close all

clear global RunData
clear global handles

global handles

scrnsize = get(0, 'screensize');
winsize = [.04*scrnsize(3) .04*scrnsize(4) floor(.9*scrnsize(3)) floor(.85* ...
    scrnsize(4))];

%create main window
handles.window = figure('units', 'pixels', 'position', winsize,...
    'name', 'ECE_527_SemesterProject_Matched Filter', 'menubar', ...
    'figure', 'numbertitle', 'off','ToolBar','figure');

set(handles.window, 'units', 'normalized');

handles.runbutton = uicontrol('Style', 'pushbutton','String', 'Run','units', 'normalized',...
    'Position', [.925 .04 .05 .05], 'Callback', 'ECE_527_SemesterProject');

handles.NUM_Targetstxt = uicontrol('Style', 'text', 'String', 'Number of Targets','units', 'normalized',...
    'Position', [.925 .19 .05 .06]);

handles.NUM_Targets = uicontrol('Style', 'edit', 'String', '2','units', 'normalized',...
    'Position', [.925 .15 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.Band_WDTtxt = uicontrol('Style', 'text', 'String', 'Chirp BW (MHz)','units', 'normalized',...
    'Position', [.925 .34 .05 .06]);

handles.Band_WDT = uicontrol('Style', 'edit', 'String', '50','units', 'normalized',...
    'Position', [.925 .30 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.Trans_Powertxt = uicontrol('Style', 'text', 'String', 'Tx Power (MW)','units', 'normalized',...
    'Position', [.925 .49 .05 .06]);

handles.Trans_Power = uicontrol('Style', 'edit', 'String', '1','units', 'normalized',...
    'Position', [.925 .45 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.Pulse_Widthtxt = uicontrol('Style', 'text', 'String', 'Pulse Width (uS)','units', 'normalized',...
    'Position', [.925 .64 .05 .06]);

handles.Pulse_Width = uicontrol('Style', 'edit', 'String', '1','units', 'normalized',...
    'Position', [.925 .60 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.Noise_Pwrtxt = uicontrol('Style', 'text', 'String', 'Noise Power (mW)','units', 'normalized',...
    'Position', [.925 .79 .05 .06]);

handles.Noise_Pwr = uicontrol('Style', 'edit', 'String', '10','units', 'normalized',...
    'Position', [.925 .75 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.NOM_Rangetxt = uicontrol('Style', 'text', 'String', 'Nominal Range(km)','units', 'normalized',...
    'Position', [.925 .94 .05 .06]);

handles.NOM_Range = uicontrol('Style', 'edit', 'String', '1','units', 'normalized',...
    'Position', [.925 .90 .05 .03], 'Callback', 'ECE_527_SemesterProject');

handles.Ant_Gaintxt = uicontrol('Style', 'text', 'String', 'Anten. Gain (dB)','units', 'normalized',...
    'Position', [.03 .19 .05 .06]);

handles.Ant_Gain = uicontrol('Style', 'edit', 'String', '20','units', 'normalized',...
    'Position', [.03 .15 .05 .03], 'Callback', 'ECE_527_SemesterProject');
handles.Smearing = uicontrol('Style', 'text', 'String', 'Smearing Coefficient','units', 'normalized',...
    'Position', [.03 .34 .05 .06]);

handles.Smearing = uicontrol('Style', 'edit', 'String', '1','units', 'normalized',...
    'Position', [.03 .30 .05 .03], 'Callback', 'radar_equation');
