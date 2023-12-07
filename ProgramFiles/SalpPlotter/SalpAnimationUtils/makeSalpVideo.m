function makeSalpVideo(salp,forcedef,ts,simResults,vidName)

    %Framerate from salp definition
    FPS = salp.visualization.FrameRate;
    dt = 1/FPS;

    %Open figure to use to make video
    newfig = figure();
    newAx = gca;
    set(newfig,'color','w');

    %Initialize video writer
    v = VideoWriter(vidName,'MPEG-4');
    v.FrameRate = FPS;
    open(v);

    %For every iteration of simulation results, plot it and add frame to
    %video
    for i = 1:size(simResults,2);

        pose = simResults(1:3,i);
        shape = simResults(4:end,i);

        plotSalpFrame(newAx,salp,pose,shape,forcedef,ts(i));

        vidFrame = getframe(newfig);
        writeVideo(v,vidFrame);

    end

    close(v);
    close(newfig);

end