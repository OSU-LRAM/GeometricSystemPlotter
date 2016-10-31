function deltritest

	fig8circ2 = @(t) (repmat([1.9;1.9],size(t))+ [cos(-pi/4) -sin(-pi/4); sin(-pi/4) cos(-pi/4)]*[.5*(1-cos(2*t)).*sign(t-pi);.5*sin(2*t)])';
	load sysplotter_data/sysf_honey_swimmer__shchf_LowRe_MaxDisplacement	

%	xy = fig8circ2(linspace(0,2*pi));
	xy = p.phi_locus_full{1}.shape;
	
	C = [(1:size(xy,1))' [2:size(xy,1) 1]' ];
	


	
	tic
	dt = DelaunayTri(xy(:,1), xy(:,2), C);
	io = dt.inOutStatus();	
	figure(101)
	clf(101)
	patch('faces',dt(io,:), 'vertices', dt.X, 'FaceColor','r');
	axis equal;
	
	Triangulation = dt.Triangulation(io,:);
	X = dt.X;
	
	
% 	i = 1;
% 	q = doubleintegral(@(a,b) ones(size(a))...
% 		,struct('type','polygon','x',X(Triangulation(i,:),1),'y',X(Triangulation(i,:),2))...
% 		,struct('method','gauss','tol',1e-6))
	
	area = tri_orient_check(Triangulation) .* arrayfun(@(i)...
		doubleintegral(@(a,b) interp2(s.grid.eval{2},s.grid.eval{1}, permute(s.height_optimized{2}, [2 1]),a,b,'cubic')...
		,struct('type','polygon','x',X(Triangulation(i,:),1),'y',X(Triangulation(i,:),2))...
		,struct('method','gauss','tol',1e-6)),(1:length(Triangulation))');
	
	a = sum(area)
	
	toc


end