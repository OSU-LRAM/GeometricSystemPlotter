%Plotter for force definition
function plotForceProfile(ax,forcedef)

    cla(ax);

    %Number of points to plot over time
    nForcePts = 1000;
    duration = forcedef.duration;

    ts = linspace(0,duration,nForcePts);
    f1 = zeros(size(ts));

    for i = 1:nForcePts
        amp = forcedef.amplitudes(ts(i));
        f1(i) = amp(1);
    end

    plot(ax,ts,f1);
    title(ax,'Force Profile for Thruster 1');
    legend(ax,'off');
    set(ax,'XLim',[0,duration]);

end