function [slines] = get_line_properties_templates(uid, style)
%GET_LINE_PROPERTIES_TEMPLATES Returns line properties for plots
%   [slines] = GET_LINE_PROPERTIES_TEMPLATES(uid, style)
%   returns line properties for plots of template attack results.
%
%   This method may be useful to get consitent line properties for each
%   compression method and parameters so that plots with same parameters
%   may use the same color, width, etc.
%
%   uid should be an integer that uniquely identifies the desired line
%   properties. This could be used to make sure that the same properties
%   are used across different figures. This method currently supports only
%   the integers from 1 to 6.
%
%   slines is a structure containing the following information:
%   - 'Color'
%   - 'LineStyle'
%   - 'LineWidth'
%   - 'Marker'
%
%   style is a string that specifies a kind of style, that may change the
%   properties returned. Currently supported styles are:
%   - 'normal'
%   - 'fancy'
%   - 'nocolor' (use for black & white images)
%   - 'cmapjet' (use for arbitrary number of lines, colormap=jet)
%
%   The slines structure can be used with make_figures_ge to use the
%   selected properties.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

%% Check and initialize parameters

%% Return line parameters
if strcmp(style, 'normal')
    slines.LineWidth = 2;
    if uid == 1
        slines.Color = 'm';
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 2
        slines.Color = 'c';
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 3        
        slines.Color = 'g';
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 4
        slines.Color = 'k';
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 5
        slines.Color = 'b';
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 6
        slines.Color = 'r';
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 7
        slines.Color = 'm';
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 8
        slines.Color = 'c';
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 9
        slines.Color = 'g';
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 10
        slines.Color = 'k';
        slines.LineStyle = '--';
        slines.Marker = 'none';  
    elseif uid == 11
        slines.Color = 'b';
        slines.LineStyle = '-.';
        slines.Marker = 'none';  
    elseif uid == 12
        slines.Color = 'r';
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 13
        slines.Color = 'y';
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    end
elseif strcmp(style, 'fancy')
    slines.LineWidth = 2;
    if uid == 1
        slines.Color = 'm';
        slines.LineStyle = '-';
        slines.Marker = 'o';
    elseif uid == 2
        slines.Color = 'c';
        slines.LineStyle = '--';
        slines.Marker = '+';
    elseif uid == 3
        slines.Color = 'g';
        slines.LineStyle = '-.';
        slines.Marker = '*';
    elseif uid == 4
        slines.Color = 'k';
        slines.LineStyle = '-';
        slines.Marker = '.';
    elseif uid == 5
        slines.Color = 'b';
        slines.LineStyle = '--';
        slines.Marker = 'x';
    elseif uid == 6
        slines.Color = 'r';
        slines.LineStyle = '-.';
        slines.Marker = 's';
    elseif uid == 7
        slines.Color = 'm';
        slines.LineStyle = '--';
        slines.Marker = '*';
    elseif uid == 8
        slines.Color = 'c';
        slines.LineStyle = '-.';
        slines.Marker = '.';
    elseif uid == 9
        slines.Color = 'g';
        slines.LineStyle = '-';
        slines.Marker = 'x';
    elseif uid == 10
        slines.Color = 'k';
        slines.LineStyle = '--';
        slines.Marker = 's';  
    elseif uid == 11
        slines.Color = 'b';
        slines.LineStyle = '-.';
        slines.Marker = '*';  
    elseif uid == 12
        slines.Color = 'r';
        slines.LineStyle = '--';
        slines.Marker = '.';
    elseif uid == 13
        slines.Color = 'y';
        slines.LineStyle = '-.';
        slines.Marker = 'x';
    end
elseif strcmp(style, 'nocolor')
    maxuid = 12;
    if uid > maxuid || uid < 1
        error('Bad uid: %d', uid);
    end
    slines.LineWidth = 2;
    cmap = colormap('Hot');
    slines.Color = cmap(1+floor((uid-1)*63/(maxuid-1)),:);
    if uid == 1
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 2
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 3        
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 4
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 5
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 6
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 7
        slines.LineStyle = '--';
        slines.Marker = 'none';
    elseif uid == 8
        slines.LineStyle = '-.';
        slines.Marker = 'none';
    elseif uid == 9
        slines.LineStyle = '-';
        slines.Marker = 'none';
    elseif uid == 10
        slines.LineStyle = '--';
        slines.Marker = 'none';  
        slines.Color = cmap(-2+floor((uid-1)*63/(maxuid-1)),:);
    elseif uid == 11
        slines.LineStyle = '-.';
        slines.Marker = 'none';  
        slines.Color = cmap(-2+floor((uid-1)*63/(maxuid-1)),:);
    elseif uid == 12
        slines.LineStyle = '-';
        slines.Marker = 'none';
        slines.Color = cmap(-2+floor((uid-1)*63/(maxuid-1)),:);
    end
elseif strcmp(style, 'cmapjet')
    maxc = 9;
    cmap = colormap('Jet');
    i = mod(uid, maxc);
    if i == 0
        i = maxc;
    end
    slines.Color = cmap(1+floor((i-1)*63/(maxc-1)),:);
    slines.LineWidth = 1 + floor((uid-1)/maxc);
    if mod(uid,3) == 0
        slines.LineStyle = '-';
    elseif mod(uid,3) == 1
        slines.LineStyle = '--';
    else
        slines.LineStyle = '-.';
    end
    slines.Marker = 'none';
else
    error('Unknown style');
end

end
