%Curve definition for 4-mode piecewise curvature
function output = curv_PiecewiseCC_4Sections(params,mode)

    %Handle call appropriately given the curvature function below
    output = curvEvalModeHandler(params,mode,@curve_fun);

end

%Define curvature along the salp spine based on shape modes
function curvature = curve_fun(s,r1,r2,r3,r4)

    %Handle any size of arclength vector
    curvature = zeros(size(s));
    %Move arclength to be in range [0,1] so I stop having a headache
    s = s + .5;

    %First section of curvature
    curvature(s<1/4) = r1;
    %Second section of curvature
    curvature(s>1/4 & s<1/2) = r2;
    %Third section of curvature
	curvature(s>1/2 & s<3/4) = r3;
    %Fourth section of curvature
    curvature(s>3/4) = r4;

end