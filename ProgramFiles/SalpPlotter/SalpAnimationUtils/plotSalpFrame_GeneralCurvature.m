%Plots salp freezeframe to a chosen axis
function plotSalpFrame_GeneralCurvature(ax,salp,pose,shape,forcedef,t)

    %Clear axis to prep for incoming drawings
    cla(ax);
    axis(ax,'equal');
    hold(ax,'on');

    base_anim_geom = generateBaseSalp_GeneralCurvature(salp,pose,shape);
    base_anim_geom.overrideAxisLimits = false;

    %If force is not defined initialize some necessary values to zero
    if nargin < 5
        t = 0;
        forcedef.amplitudes = @(t) zeros(1,salp.geometry.thrusters.nThrusters);
        forcedef.maxThrust = 1;
        base_anim_geom.overrideAxisLimits = true;
    end
        


    plotGrid(ax,base_anim_geom);
    plotCurrentSalpGeometry_GeneralCurvature(ax,base_anim_geom,salp,pose,forcedef,t);

    formatSalpPlot(ax,salp,base_anim_geom);

    drawnow;

end