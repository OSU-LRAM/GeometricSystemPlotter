%Plots salp freezeframe to a chosen axis
function plotSalpFrame_NLinkChain(ax,salp,pose,shape,forcedef,t)

    %Clear axis to prep for incoming drawings
    cla(ax);
    axis(ax,'equal');
    hold(ax,'on');

    %Get structcure holding salp configuration information
    base_anim_geom = generateBaseSalp_NLink(salp,pose,shape);
    base_anim_geom.overrideAxisLimits = false;

    %If force not defined initialize some necessary force information to
    %zero
    if nargin < 5
        t = 0;
        forcedef.amplitudes = @(t) zeros(1,salp.geometry.thrusters.nThrusters);
        forcedef.maxThrust = 1;
        base_anim_geom.overrideAxisLimits = true;
    end
        
    %Put the grey grid in the background
    plotGrid(ax,base_anim_geom);
    %Plot the salp at current configuration
    plotCurrentSalpGeometry_NLinkChain(ax,base_anim_geom,salp,pose,forcedef,t);

    formatSalpPlot(ax,salp,base_anim_geom);

    drawnow;

end