function B = fatbackbone_from_curvature_bases(kappa_basis_input,r,L,width)

	[h, J] = backbone_from_curvature_bases(kappa_basis_input,r,L);
	
	
	B = fatbackbone(h,L*[-.5 .5],width);


end