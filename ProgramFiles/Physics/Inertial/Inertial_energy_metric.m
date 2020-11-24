function Mp = Inertial_energy_metric(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

% Identify what kind of system is being calculated, and use this to specify how
% the local connection should be generated
switch geometry.type
    
    case {'curvature basis','curvature bases','general curvature'}
        physics_function = @Inertial_metric_continuous;
        
    case {'n-link chain','branched chain'}
        physics_function = @Inertial_metric_discrete;
        
end

% Call the physics function identified for the system
Mp = physics_function(geometry,physics,shapeparams);



end