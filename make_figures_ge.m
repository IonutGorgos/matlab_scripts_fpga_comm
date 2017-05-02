function make_figures_ge(G, xvec, ...
                         rpath, rprefix, ...
                         title_ge, slegend, font_size, ...
                         slines, options, yrange, ...
                         bounds, xlstr, ylstr)
%MAKE_FIGURES_GE creates figures for guessing entropy data
%   MAKE_FIGURES_GE(G, xvec, ...
%                        rpath, rprefix, ...
%                        title_ge, slegend, font_size, ...
%                        slines, options, yrange, bounds, xlstr, ylstr)
%   creates figures for guessing entropy data.
%
%   G must be a matrix of size nr_exp x nx, where nr_exp
%   specifies the number of experiments and nx is the
%   number of different values (to be ploted along the x axis) for which
%   the guessing entropy was computed at each experiment.
%
%   As an example, the success_info structure returned by
%   get_success_info_like can be used to produce guessing entropy values
%   using the function get_ge_from_success_info.
%
%   xvec is a vector containing the data for the x-axis in figures.
%
%   rpath is a string containing the path where the figures should be
%   saved.
%
%   rprefix is a string that will be used to prefix all saved figures.
%
%   title_result specifies an overall title for the figures. Pass []
%   (empty) to ignore.
%
%   slegend should be a cell of strings containing the legend to be shown
%   for each experiment. Pass [] to omit legends.
%
%   font_size specifies the font size of the text. Pass [] to use the
%   default (24).
%
%   slines should be a cell of structures of length nr_exp, containing line
%   properties to be used with each experiment. Each structure in slines
%   should contain the following fields:
%   - 'Color'
%   - 'LineStyle'
%   - 'LineWidth'
%   - 'Marker'
%   See "doc plot" for details on these parameters.
%
%   options is a string that can specify options for the plot.
%   Pass one or more of the following values if desired:
%       'y': use the default Matlab ylim values for the given data (to
%       maximize the space over which the data is shown). If this is not
%       given then the range specified by yrange (see below) is used.
%       'g': use grid on. Default is off.
%       'L': use large figure size (useful with 'Publish' mode).
%       'p': print PDF for publication.
%       'i': Save a PNG image.
%       'b': use the 'boundedline' script to produce bounded lines (useful
%       to show confidence regions, etc.).
%       't': use 'XTick' and 'xlim' from given xvec.
%       'l': use logarithmic scale for the Y-axis.
%       'n': use linear scale for the X-axis.
%       'f': do not append any string to given file prefix when saving fig.
%       's': put the legend in the south-east corner rather than north-east
%       '2': add a "2^" in front of each YTickLabel. Useful when ploting
%       data that is in log2. Plot Y axis in log2 rather than log10.
%
%   yrange is a 2-element vector that specifies the ylim values for the
%   figures.
%
%   bounds should be a matrix of size nr_exp x nr_test_groups, specifying
%   the bounds (confidence region, standard dev, standard error, etc) for
%   each of the lines given to plot. This is an optional parameter. If
%   given, together with option 'b', then the 'boundedline' script is
%   used to create a plot with bounded shaded lines. Colors in slines are
%   ignored for now. See this link:
%   http://www.mathworks.com/matlabcentral/fileexchange/
%   27485-boundedline-line-plots-with-shaded-errorconfidence-intervals
%
%   xlstr, if given, will be used as xlabel (instead of default '$n_a$').
%
%   ylstr, if given, will be used as ylabel (instead of default 'guessing_entropy').

%% Initialize and check stuff
nr_exp = size(G, 1);
if isempty(font_size)
    font_size = 12;
end
nx = length(xvec);
if size(G, 2) ~= nx
    error('Incompatible nx');
end
font_small = font_size - 4;
if strfind(options, 'L')
    fig_rect = [1 1 1024 768];
else
    fig_rect = [1 1 640 480];
end
if nargin < 11
    bounds = [];
end
if nargin < 12
    xlstr = '$n_a$ (log axis)';
end
if nargin < 13
    ylstr = 'Guessing entropy (bits)';
end

%% Plot guessing entropy data
h = figure('Position', fig_rect);

if ~isempty(strfind(options, 'p'))
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [12 10]);
    set(gcf, 'PaperPosition', [0.2 0.2 11.7 9.7]);
end

if ~isempty(strfind(options, 'b')) && ~isempty(bounds)
    [hl, ~] = boundedline(xvec, G, permute(bounds, [2, 3, 1]));
    set(gca, 'XScale', 'log');
    for i=1:nr_exp
        set(hl(i), 'LineStyle', slines{i}.LineStyle);
        set(hl(i), 'LineWidth', slines{i}.LineWidth);
        set(hl(i), 'Marker', slines{i}.Marker);
    end
else
    for i=1:nr_exp
        if (~isempty(strfind(options, 'l')) && ~isempty(strfind(options, '2')))
            semilogx(xvec, log2(G(i,:)), ...
            'Color', slines{i}.Color, ...
            'LineStyle', slines{i}.LineStyle, ...
            'LineWidth', slines{i}.LineWidth, ...
            'Marker', slines{i}.Marker); 
        else
            semilogx(xvec, G(i,:), ...
                'Color', slines{i}.Color, ...
                'LineStyle', slines{i}.LineStyle, ...
                'LineWidth', slines{i}.LineWidth, ...
                'Marker', slines{i}.Marker);   
        end
        hold on;
    end
    hold off;
end

if isempty(strfind(options, 'y'))
    ylim(yrange);
end

set(gca,'FontSize', font_size);
title_str = sprintf('%s', title_ge);
set(0, 'DefaulttextInterpreter', 'latex');

if ~isempty(xlstr)
    lh = xlabel(xlstr);
    set(lh, 'FontSize', font_size);
end

if ~isempty(title_ge)
    title(title_str);
end
set(0, 'DefaulttextInterpreter', 'none');

if ~isempty(ylstr)
    lh = ylabel(ylstr);
    set(lh, 'FontSize', font_size);
end

if ~isempty(slegend)
    lh = legend(slegend);
    if ~isempty(strfind(options, 's'))
        set(lh, 'Location', 'SouthEast');
    end
    set(lh, 'interpreter', 'none');
    set(lh, 'FontSize', font_small);
end

if ~isempty(strfind(options, 'g'))
    grid on;
end

if isempty(strfind(options, 'p'))
    orient landscape;
end

if ~isempty(strfind(options, 't'))
    set(gca, 'XTick', xvec);
    set(gca, 'xlim', [xvec(1), xvec(end)]);
end

if ~isempty(strfind(options, 'l'))
    if isempty(strfind(options, '2'))
        set(gca,'YScale', 'log');
    else
        yt = get(gca, 'YTick');
        nr_y_ticks = length(yt);
        ytl = cell(nr_y_ticks, 1);
        for j=1:nr_y_ticks
            ytl{j} = ['2^' num2str(yt(j))];
        end
        set(gca, 'YTickLabel', ytl)
    end
end

if ~isempty(strfind(options, '2'))
    
end

if ~isempty(strfind(options, 'n'))
    set(gca,'XScale', 'linear');
end

if ~isempty(strfind(options, 'f'))
    fname = sprintf('%s%s', rpath, rprefix);
else
    fname = sprintf('%s%sguess_entropy', rpath, rprefix);
end

if ~isempty(strfind(options, 'p'))
    print(h, '-dpdf', fname);
end

if ~isempty(strfind(options, 'i'))
    print(h, '-dpng', fname);
end

end
