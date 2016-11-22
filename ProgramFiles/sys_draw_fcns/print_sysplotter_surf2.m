function print_sysplotter_surf2(f,destination)
%wrapper for the painter_layer function, with arguments for the surface
%plots used in sysplotter

%Set the figure to print with no background
ihc = get(f,'InvertHardcopy');
c = get(f,'Color');
set(f,'InvertHardCopy','off','Color','none')

%get the handle of the target axes
ax = findobj(f,'Type','axes');


%define layers and send for printing
try
painter_layer2(f,...
	{ax,findobj(ax,'Type','surface'),findobj(ax,'Type','patch')},...
	{'vector','raster','vector'},...
	destination,...
	{'-noui'});
catch
	%
end

%restore the figure properties
set(f,'InvertHardCopy',ihc,'Color',c)