function sysprint2(f,destination)
%Print a figure created by sysplotter

%Get the userdata for the figure
udata = get(f,'UserData');

switch udata.plottype
	
	case {'vfield','hfun-contour','BVI','disp','dbeta'}
		
		flatprint(f,destination,{'-noui','-loose'});
		
	case {'hfun-surface','beta'}
		
		print_sysplotter_surf2(f,destination);
		
	otherwise
		
		error('Unknown figure type')
		
end