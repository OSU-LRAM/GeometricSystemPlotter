%Decides which force connection to use for salp body velocity calculation
%depending on system physics regime
function [g_circ,shapeVel] = salpForceToVelocity(geometry,physics,shape,forces)

    switch physics.systemType
        case 'drag'
            [g_circ,shapeVel] = salpForceToVelocity_DragDominated(geometry,physics,shape,forces);
        case 'inertia'
            error('Inertia-dominated swimming not yet implemented for salps');
    end

end