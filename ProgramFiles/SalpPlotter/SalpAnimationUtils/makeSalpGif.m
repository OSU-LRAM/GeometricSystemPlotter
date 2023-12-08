function makeSalpGif(salp,forcedef,ts,simResults,gifName)

    %Get framerate from salp definition
    FPS = salp.visualization.FrameRate;
    dt = 1/FPS;

    %Open window to use to make gif
    newfig = figure();
    newAx = gca;
    set(newfig,'color','w');

    %For every sim time, plot salp at given pose
    for i = 1:size(simResults,2);

        pose = simResults(1:3,i);
        shape = simResults(4:end,i);

        plotSalpFrame(newAx,salp,pose,shape,forcedef,ts(i));

        if i == 1
            gif(gifName,'DelayTime',dt,'overwrite',true);
        else
            gif;
        end
    end

    close(newfig);

end