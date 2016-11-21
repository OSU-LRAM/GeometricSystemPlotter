function tri_sign = tri_orient_check(tri)
% Check the orientation of a set of triangles with respect to the ordering
% of their nodes (Delaunay triangles are always positively oriented
% geometrically. This test checks if this ordering coincides with the order
% of the points in the polygon.)

	tri_sign = (tri(:,1) < tri(:,2)) & (tri(:,2) < tri(:,3))... % Get all positive orientations
		| (tri(:,2) < tri(:,3)) & (tri(:,3) < tri(:,1)) ...
		| (tri(:,3) < tri(:,1)) & (tri(:,1) < tri(:,2));

	tri_sign = -1+2*tri_sign; % Make all zeros -1

end