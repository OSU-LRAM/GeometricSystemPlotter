function footprint = fatbackbone(backbone_fun,endpoints,width)
% Make a polygon from a backbone curve

	% get the positions and orientations of points along the core
	corepoints = linspace(endpoints(1),endpoints(2),50);

	core = backbone_fun(corepoints);
	
	% make primitives for the end
	a = .3*range(endpoints)/3/2;  % to match end of three-link system
	b = width/2;

	poscaptheta = linspace(pi/2,-pi/2,30)';
	poscapprim = [a*cos(poscaptheta) b*sin(poscaptheta)];
	
	negcaptheta = linspace(3*pi/2,pi/2,30)';
	negcapprim = [a*cos(negcaptheta) b*sin(negcaptheta)];
	
	
	
	% find the location of points along the perimeter
	posedge = zeros(length(corepoints),2);
	negedge = posedge;
	

	
	for i = 1:length(corepoints)
		
		x = core(1,i);
		y = core(2,i);
		th = core(3,i);
		
		
		transform = [cos(th) -sin(th) x;
			sin(th) cos(th) y;
			0 0 1];
		
		edgepoints = transform*[0 0;width/2 -width/2;1 1];
		
		posedge(i,:) = edgepoints(1:2,1)';
		negedge(end+1-i,:) = edgepoints(1:2,2)';
		
		if i == 1
			negcap = (transform*[negcapprim';ones(1,length(negcaptheta))])';
			negcap(:,3) = [];
		end
		
		if i == length(corepoints)
			poscap = (transform*[poscapprim';ones(1,length(poscaptheta))])';
			poscap(:,3) = [];
		end

		
	end	
	
	

	footprint = [posedge;poscap;negedge;negcap];
	

end