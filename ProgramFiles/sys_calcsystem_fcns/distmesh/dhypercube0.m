function [d d_component closest_face] = dhypercube0(p,negcorner,poscorner)

	% Replicate the corner values to be the same size as the point array
	negcorner = repmat(negcorner,[size(p,1), 1]);
	poscorner = repmat(poscorner,[size(p,1), 1]);
	
	% Get the signed component-wise displacements from the points to the
	% corners. These correspond to distances from the (n_dim-1)-faces of
	% the hypercube
	d1 =  negcorner-p;
	d2 = -poscorner+p;
	
	% save out the distances for use separate use
	d_component = -min(-d1,-d2);
	
	% Mark which face this distance is to
	closest_face = ones(size(d_component));
	closest_face(d_component == d1) = -1;
	
	
	% Find the minimum-length components of the displacements
	d1min = -min(-d1,[],2);
	d2min = -min(-d2,[],2);
	
	% Get the negative minimum distance overall
	d = -min(-d1min,-d2min);
	
	
	% For points outside the hypercube, use the L2-norm for over all outside
	% components
	
	% Indices of points outside the hypercube
	oI = any(d1>0,2) | any(d2>0,2);
	
	% Signed componentwise distances of those points
	d1o = d1(oI,:);
	d2o = d2(oI,:);
	
	% Mask out the components for which the point is "inside" the hypercube
	d1o = d1o.*(d1o>0);
	d2o = d2o.*(d2o>0);
	
	% L2 norm applied to those points
	d(oI,:) = sqrt(sum( (d1o.^2) + (d2o.^2) ,2) );
	
end