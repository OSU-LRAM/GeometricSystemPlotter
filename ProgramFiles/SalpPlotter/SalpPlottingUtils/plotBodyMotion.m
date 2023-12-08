%Plotter for body motion over time from simulation data
function plotBodyMotion(ax,simResults)

    cla(ax);

    plot(ax,simResults(1,:),simResults(2,:));

    minX = min(simResults(1,:))-1;
    maxX = max(simResults(1,:))+1;
    minY = min(simResults(2,:))-1;
    maxY = max(simResults(2,:))+1;

    axis(ax,[minX,maxX,minY,maxY]);
    title(ax,'Body Motion Results');
    axis(ax,'equal');

end