% metric value checker for debugging.

% Inputs
clearvars -except sysplotter_inputpath

snake_sys = 'serpenoid';    % Which backbone to use?
% snake_sys = 'piecewise';

lambda_func = 'quad';       % Which lambda input function?
% lambda_func = 'abs';

if strcmp(lambda_func,'quad')
    b = round(linspace(0,0.101,11),3); % quadratic value b range
%     b = b(11:15);
elseif strcmp(lambda_func,'abs')
    b = round(linspace(0,0.318,21),3); % absolute value b range
    %         b = b(1);
end

plot_idx = [1,11:15];
tic
for i = 1:numel(b)
    
        if strcmp(snake_sys,'serpenoid') % serpenoid system
            if strcmp(lambda_func,'quad') % quadratic lambda function
                current_system = ['sysf_serpenoid_extendable_',num2str(b(i)*1000)];
                curvfilename = ['serpenoid_stretch_b_',num2str(b(i)*1000)];
            elseif strcmp(lambda_func,'abs') % absolute lambda function
                current_system = ['sysf_serpenoid_extendable_abs_',num2str(b(i)*1000)];
                curvfilename = ['serpenoid_stretch_abs_b_',num2str(b(i)*1000)];
            end
        elseif strcmp(snake_sys,'piecewise') % piecewise system
            if strcmp(lambda_func,'quad')
                current_system = ['sysf_piecewise_const_stretch_b_',num2str(b(i)*1000)];
                curvfilename = ['piecewise_const_stretch_b_',num2str(b(i)*1000)];
            elseif strcmp(lambda_func,'abs')
                current_system = ['sysf_piecewise_const_stretch_abs_b_',num2str(b(i)*1000)];
                curvfilename = ['piecewise_const_stretch_abs_b_',num2str(b(i)*1000)];
            end
        end
    disp(['    system: ',current_system])
    
    load([current_system,'.mat'])
    
    curv_def = ['curv_',curvfilename,'.m'];
%     metric_eval = s.density.metric_eval;
    eval_density = s.density.eval;
    grid_range = s.grid_range;
    
    x_eval = linspace(grid_range(1),grid_range(2),eval_density(1));
    y_eval = linspace(grid_range(3),grid_range(4),eval_density(2));
    
    [x,y] = ndgrid(x_eval,y_eval);
    metric_det = zeros(eval_density(1),eval_density(2));
    
    disp(['    Gathering Metric Determinants for system ',current_system])
    for j = 1:eval_density(1)
        for k = 1:eval_density(2)
            
            metric_det(j,k,i) = det(s.metric(x_eval(j),y_eval(k)));
            
        end
    end
    disp('    Done.')
    disp('------------')
end
disp(['    Finished evaluating ',snake_sys,' ',lambda_func,' in ',toc/60,' minutes.'])
disp('Plotting desired data.')

for i = 1:numel(plot_idx)

    z = metric_det(:,:,plot_idx(i));
    figure(i); surf(x,y,z)
    title(['Metric determinant, ',snake_sys,' ',lambda_func,' b ',num2str(b(plot_idx(i)))]);

end
