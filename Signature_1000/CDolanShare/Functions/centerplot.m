function [ output_args ] = centerplot( width, height, unit )
%CENTERPLOT Re-sizes and centers plot based on specified size
%   CENTERPLOT centers the current figure
%   
%   CENTERPLOT(width,height) centers the current figure and resizes the
%   figure using the same units that are used by the figure
%
%   CENTERPLOT(width,height,unit) centers the current figure and resizes
%   the figure using the specified units (only 'inches', 'centimeters', or
%   'pixels' can be used).
%
%   A.P. Garcia

%%

% Maintains same size if width and height are not identified
if nargin < 2
    pos = get(gcf,'Position');
    width = pos(3);
    height = pos(4);
end

% Maintains same units if units are not identified
if nargin < 3
    unit = get(gcf,'Units');
end

set(gcf,'Units',unit)

% Find current figure unit, convert to pixels
figunit = get(gcf,'Units');
if isequal(figunit,'inches')
    figscale = get(0,'ScreenPixelsPerInch');
elseif isequal(figunit,'pixels')
    figscale = 1;
elseif isequal(figunit,'centimeters')
    figscale = get(0,'ScreenPixelsPerInch')/2.54;
end

% Find scaling factor for new unit, set to pixels
if isequal(unit,'inches')
    scale = get(0,'ScreenPixelsPerInch');
elseif isequal(unit,'pixels')
    scale = 1;
elseif isequal(unit,'centimeters')
    scale = get(0,'ScreenPixelsPerInch')/2.54;
end

width = width*scale;
height = height*scale;

% Gets screen size
Screen = get(0,'ScreenSize');
NewPos = [(Screen(3)-width)/2,(Screen(4)-height)/2,width,height]/figscale;
    
% Sets new size of figure for screen as well as for exporting
set(gcf,'PaperPosition',NewPos)
set(gcf,'Position',get(gcf,'PaperPosition'))

end

