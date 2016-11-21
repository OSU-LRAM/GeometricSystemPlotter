function hiressurfprep(f)

load sys_draw_fcns/BlackWhiteRedColormap
set(f, 'colormap',bwr)

ax = findobj(f,'type','axes');
axis(ax,'vis3d')
set(ax,'DataAspectRatio',[30,30,1])
view(ax,16,50)
nicebox(ax,'redraw')

gait = findobj(ax,'type','patch');
set(gait,'linewidth',6)

set(ax,'XtickLabel',[])
set(ax,'YtickLabel',[])
set(ax,'ZtickLabel',[])
set(ax,'TickLength',[0 0])
xlabel(ax,'')
ylabel(ax,'')

sysprint(f,'~/Temp/tempsurface.pdf')

axis(ax,'off')
nicebox(ax,'off')

%%%
%Upsample surface
upsamplefactor=3;
s = findobj(ax,'type','surface');

xdata = get(s,'xdata');
ydata = get(s,'ydata');
zdata = get(s,'zdata');

xrange = [min(xdata(:)),max(xdata(:))];
yrange = [min(ydata(:)),max(ydata(:))];
zrange = [min(zdata(:)),max(zdata(:))];

[xnew,ynew] = meshgrid(linspace(xrange(1),xrange(2),upsamplefactor*size(xdata,1)),...
	linspace(yrange(1),yrange(2),upsamplefactor*size(ydata,2)));

znew = interp2(ydata,xdata,zdata,ynew,xnew,'cubic');




set(s,'xdata',xnew,'ydata',ynew,'zdata',znew)


set(f,'PaperSize',[8.5 11]*8)
set(f,'Paperposition',get(gcf,'Paperposition')*8)
set(f,'color',[0 1 0])
shading(ax,'interp')
print(f,'-dpng', '~/Temp/greenscreensurface')


% xgrid = [xnew([(1:3:end) end],:)',xnew(:,[(1:3:end) end])];
% xgrid = [xgrid;NaN(1,size(xgrid,2))];
% ygrid = [ynew([(1:3:end) end],:)',ynew(:,[(1:3:end) end])];
% ygrid = [ygrid;NaN(1,size(ygrid,2))];
% zgrid = [znew([(1:3:end) end],:)',znew(:,[(1:3:end) end])];
% zgrid = [zgrid;NaN(1,size(zgrid,2))];
% 
% f2 = figure;
% set(f2,'Color','None','InvertHardCopy','off')
% ax2 = axes('parent',f2);
% axis(ax2,'vis3d')
% set(ax2,'DataAspectRatio',[30,30,1])
% view(ax2,16,50)
% axis(ax2,'off')
% 
% line(xgrid(:),ygrid(:),zgrid(:),'Color',[1 1 1]*0.3,'LineWidth',0.15);

%set(s,'visible','off')
%set(gait,'visible','off')