%  new data plotter. Only going for contour or metric plots, just need
%  something simpler than the old one.

clearvars -except sysplotter_inputpath


%%% system parameters
% snake_sys_array = {'serpenoid','piecewise'};  % Which backbone to use?
snake_sys_array = {'serpenoid'};

% lambda_func_array = {'quad','abs'};         % Which lambda input function?
lambda_func_array = {'quad'};

% length_sys_array = {'original','Lrms'};       % Lrms length or changeable?
length_sys_array = {'original'};


%%% plot display parameters
% metric_stretch_array = {'no','metric'};
metric_stretch_array = {'no'};

% plot_gait_array = {'with gait','no gait'};
plot_gait_array = {'with gait'};

save_plot = 0;

extra_text = '';
% extra_text = '_lg_grid';

num_b = 21;

% b = round(linspace(0,0.101,21),3);
% b = 0:0.01:0.1;
% b = 0.07;

fig = figure(3);
tic
for id1 = 1:numel(snake_sys_array)
    snake_sys = snake_sys_array{id1};
    
    for id2 = 1:numel(lambda_func_array)
        lambda_func = lambda_func_array{id2};
        
        for id3 = 1:numel(length_sys_array)
            length_sys = length_sys_array{id3};
            
            for id4 = 1:numel(metric_stretch_array)
                metric_stretch = metric_stretch_array{id4};
                
                for id5 = 1:numel(plot_gait_array)
                    plot_gait = plot_gait_array{id5};
                    
                    if strcmp(lambda_func,'quad')
                        b = round(linspace(0,0.101,num_b),3); % quadratic value b range
                        %         b = b(1);
                    elseif strcmp(lambda_func,'abs')
                        b = round(linspace(0,0.318,num_b),3); % absolute value b range
                        %         b = b(1);
                    end
                    
                    for i = 1:numel(b)
                        if strcmp(length_sys,'original')
                            if strcmp(snake_sys,'serpenoid')
                                if strcmp(lambda_func,'quad')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_serpenoid_extendable_',num2str(b(i)*1000),'__shchf_serpenoid_seed_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_serpenoid_extendable_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                elseif strcmp(lambda_func,'abs')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_serpenoid_extendable_abs_',num2str(b(i)*1000),'__shchf_serpenoid_seed_abs_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_serpenoid_extendable_abs_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                end
                            elseif strcmp(snake_sys,'piecewise')
                                if strcmp(lambda_func,'quad')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_piecewise_const_stretch_b_',num2str(b(i)*1000),'__shchf_piecewise_seed_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_piecewise_const_stretch_b_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                elseif strcmp(lambda_func,'abs')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_piecewise_const_stretch_abs_b_',num2str(b(i)*1000),'__shchf_piecewise_seed_abs_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_piecewise_const_stretch_abs_b_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                end
                            end
                        elseif strcmp(length_sys,'Lrms')
                            if strcmp(snake_sys,'serpenoid')
                                if strcmp(lambda_func,'quad')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_serpenoid_quad_Lrms_',num2str(b(i)*1000),'__shchf_serpenoid_seed_Lrms_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_serpenoid_quad_Lrms_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                elseif strcmp(lambda_func,'abs')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_serpenoid_abs_Lrms_',num2str(b(i)*1000),'__shchf_serpenoid_seed_abs_Lrms_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_serpenoid_abs_Lrms_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                end
                            elseif strcmp(snake_sys,'piecewise')
                                if strcmp(lambda_func,'quad')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_piecewise_quad_Lrms_',num2str(b(i)*1000),'__shchf_piecewise_seed_Lrms_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_piecewise_quad_Lrms_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                elseif strcmp(lambda_func,'abs')
                                    if strcmp(plot_gait,'with gait')
                                        system_file = ['sysf_piecewise_abs_Lrms_',num2str(b(i)*1000),'__shchf_piecewise_seed_abs_Lrms_b_',num2str(b(i)*1000),'_optimal.mat'];
                                    elseif strcmp(plot_gait,'no gait')
                                        system_file = ['sysf_piecewise_abs_Lrms_',num2str(b(i)*1000),'__null.mat'];
                                    end
                                end
                            end
                        end
                        
                        load(system_file)
                        
                        % Prep plot stuff
                        
                        if strcmp(metric_stretch,'no') || s.no_flatten
                            % get the contour lines
                            g1s = s.grid.eval{1,1}; % grid lines
                            g2s = s.grid.eval{2,1};
                            
                            if strcmp(plot_gait,'with gait')
                                t = p.time{1,1}{1,1}; %#ok<*UNRCH>
                                alpha = p.phi_def{1,1}{1,1}(t);
                            end
                            
                        elseif strcmp(metric_stretch,'metric')
                            
                            g = s.grid.eval;
                            H = cat(1,s.DA,s.DA_optimized); % currently ignores singularities
                            
                            % Get the value by which to scale the constraint curvature function
                            ascale = arrayfun(@(x,y) 1/det(s.convert.jacobian(x,y)),g{:});
                            
                            % Apply the jacobian to the vectors
                            H = cellfun(@(x) x.*ascale,H,'UniformOutput',false);
                            
                            % Convert the grid points to their new locations
                            [g{:}] = s.convert.old_to_new_points(g{:});
                            grid(:) = g;
                            
                            edgeres = 30;
                            
                            oldx_edge = [s.grid_range(1)*ones(edgeres,1);linspace(s.grid_range(1),s.grid_range(2),edgeres)';...
                                s.grid_range(2)*ones(edgeres,1);linspace(s.grid_range(2),s.grid_range(1),edgeres)'];
                            oldy_edge = [linspace(s.grid_range(3),s.grid_range(4),edgeres)';s.grid_range(4)*ones(edgeres,1);...
                                linspace(s.grid_range(4),s.grid_range(3),edgeres)';s.grid_range(3)*ones(edgeres,1)];
                            
                            [x_edge,y_edge] = s.convert.old_to_new_points(oldx_edge,oldy_edge);
                            
                            if strcmp(plot_gait,'with gait')
                                t = p.time{1,1}{1,1}; %#ok<*UNRCH>
                                alpha = p.phi_def{1,1}{1,1}(t);
                                
                                oldx_alpha = alpha(:,1);
                                oldy_alpha = alpha(:,2);
                                
                                [alpha(:,1),alpha(:,2)] = s.convert.old_to_new_points(oldx_alpha,oldy_alpha);
                                
                            end
                            
                            g1s = grid{1};
                            g2s = grid{2};
                            
                        end
                        z = s.DA_optimized{1,1}; % actual contour height data
                        
                        
                        mn = min(min(z));
                        mx = max(max(z));
                        lns = linspace(mn,mx,9);
                        
                        clrmap = color_Red_new(numel(lns));
                        
                        % do the plot stuff and save
                        
                        clf
                        
                        title([snake_sys,' ',lambda_func,' ',length_sys,' ',metric_stretch,' stretch ',plot_gait,' b = ',num2str(b(i))]);
                        hold on;
                        
                        
                        for j = 1:numel(lns)
                            contour(g1s,g2s,z,[lns(j) lns(j)],'Color',clrmap.colormap(j,:),'LineWidth',2);
                        end
                        contour(g1s,g2s,z,[0 0],'Color',[0 0 1],'LineWidth',2);
                        
                        if strcmp(metric_stretch,'metric')
                            l_edge = line('Parent',gca,'Xdata',x_edge,'YData',y_edge,'Color','k','LineWidth',1);
                        end
                        
                        % plot specified values
                        
                        %     % min
                        %     [min_col,Ic] = min(z);
                        %     [min_row,Ir] = min(min_col);
                        %     labelstr = sprintf('%.4f', z(Ir,Ic(Ir)));
                        % %     plot(g1s(Ir,1),g2s(1,Ic(Ir)),'k+')
                        %     text(g1s(Ir,1), g2s(1,Ic(Ir)), labelstr);
                        %
                        %     % max
                        %     [max_col,Ic] = max(z);
                        %     [max_row,Ir] = max(max_col);
                        %     max_id = [Ir Ic(Ir)];
                        %     labelstr = sprintf('%.4f', z(Ir,Ic(Ir)));
                        % %     plot(g1s(Ir,1),g2s(1,Ic(Ir)),'k+')
                        %     text(g1s(Ir,1), g2s(1,Ic(Ir)), labelstr);
                        
                        % center
%                         labelstr = sprintf('%.4f', z(11,11));
%                         plot(g1s(11,1),g2s(1,11),'k+')
%                         text(g1s(11,1), g2s(1,11),labelstr);
                        
                        if strcmp(plot_gait,'with gait')
                            plot(alpha(:,1),alpha(:,2),'Color','r','LineWidth',4)
                        end
                        axis equal
                        box on
%                         axis([-16 16 -16 16])
                        pause(0.5)
                        if save_plot % save .png and .fig files for later use
                            saveas(fig,[snake_sys,' ',lambda_func,' ',length_sys,' ',metric_stretch,' stretch ',plot_gait,' b ',num2str(b(i)*1000),extra_text,'.png']);
                            saveas(fig,[snake_sys,' ',lambda_func,' ',length_sys,' ',metric_stretch,' stretch ',plot_gait,' b ',num2str(b(i)*1000),extra_text,'.fig']);
                        end
                    end
                end
            end
        end
    end
end

clearvars -except sysplotter_inputpath

toc


