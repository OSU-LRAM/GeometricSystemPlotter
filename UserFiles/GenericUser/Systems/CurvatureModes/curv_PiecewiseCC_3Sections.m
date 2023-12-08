%Curve definition for 3 mode piecewise continuous salp
function output = curv_PiecewiseCC_3Sections(params,mode)

    %Handle calls to this function appropriately given the curve definition
    output = curvEvalModeHandler(params,mode,@curve_fun);

end

%Returns curvature at desired arclength given the shape modes
function curvature = curve_fun(s,r1,r2,r3)

    curvature = zeros(size(s));
    %Work in range [0,1] because it's easier to think about
    s = s + .5;

    curvature(s<1/3) = r1;
    curvature(s>1/3 & s<2/3) = r2;
	curvature(s>2/3) = r3;

end