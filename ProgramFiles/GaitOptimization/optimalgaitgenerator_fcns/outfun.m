function stop=outfun(y,optimValues,state,stretch,s,handles)
%%%%%%%%% 
%
%This function plots the current state of the gait on the sysplotter GUI
%after every iteration of the optimizer
%
%%%%%%%%% 

n=100;
dimension=length(y(1,:));

% % The if else statement below deletes gaits 2 iterations after they have been plotted
% if optimValues.iteration>2
%     children=get(gca,'children');
%     delete(children(6:10));
% else
% end

for thisAxes = [1:numel(handles.plot_thumbnails.Children)]
    
    axes(handles.plot_thumbnails.Children(thisAxes));

    % The if else statement below fades the gait plotted during the previous iteration
    if optimValues.iteration>1
        children=get(gca,'children');
        for idx = 1:numel(children)

            if iscell(children(idx).UserData) && strcmp(children(idx).UserData{1},'OptimizeTracer')
                children(idx).UserData = {'OptimizeTracer', children(idx).UserData{2}-1};

                if children(idx).UserData{2} == 0

                    delete(children(idx));

                else

                    children(idx).Color=[0.5 0.5 0.5];
                    children(idx).LineWidth=4;
                end
            end

        end
    %     children(1).Color=[0.5 0.5 0.5];
    %     children(2).Color=[0.5 0.5 0.5];
    %     children(3).Color=[0.5 0.5 0.5];
    %     children(4).Color=[0.5 0.5 0.5];
    %     children(5).Color=[0.5 0.5 0.5];
    % 
    %     children(1).LineWidth=4;
    %     children(2).LineWidth=4;
    %     children(3).LineWidth=4;
    %     children(4).LineWidth=4;
    %     children(5).LineWidth=4;
    else
    end

    % The if else statement below plots the gait after every iteration
    if optimValues.iteration>0
        y1 = path_from_fourier(y,n,dimension);
        hold on
        if stretch
            stretchnames = {'stretch','surface'};
            stretchname = stretchnames{stretch};

            [x_temp,y_temp,z_temp] = s.convert.(stretchname).old_to_new_points(y1(:,1),y1(:,2));
        else
            x_temp = y1(:,1);
            y_temp = y1(:,2);
            z_temp = zeros(size(y1(:,1)));
            if size(y1,2) > 2
                z_temp = y1(:,3);
            end
        end
        handle1=line('XData',x_temp,'YData',y_temp,'ZData',z_temp,'color','k','linewidth',3,'UserData',{'OptimizeTracer',2}); %#ok<NASGU>
        %plot_dir_arrows(y1(:,1),y1(:,2),2,'Color',[0 0 0],'LineWidth',3);
    else
    end

end

pause(0.05)
stop=false;
end