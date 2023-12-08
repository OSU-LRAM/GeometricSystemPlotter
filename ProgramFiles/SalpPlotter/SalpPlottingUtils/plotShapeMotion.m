%Plotter for salp shapes over time from simulation data
function plotShapeMotion(ax,ts,simResults)

    cla(ax);
    hold(ax,'on');

    nShapes = size(simResults,1)-3;

    for i = 1:nShapes
        plot(ax,ts,simResults(i+3,:),'DisplayName',['Shape ',num2str(i)]);
    end

    legend(ax);