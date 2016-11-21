function B = fatbackbone_from_general_curvature(curvdef,cparams,L,width)

	[h] = backbone_from_general_curvature(curvdef,cparams,L);
	
	
	B = fatbackbone(h,L*[-.5 .5],width);


end