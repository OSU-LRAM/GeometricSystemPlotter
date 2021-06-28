function [cBVI_to, cBVI_opt_to] = add_tlb_plot(s,p,is_opt,is_square, ignore_bounds)
    % sanitation
    if ~ishandle(1)
        error('Open trajectory graph before adding third-order effects');
    elseif ~exist('is_opt','var')
        is_opt = false;
    elseif ~exist('is_square', 'var')
        is_square = false;
    elseif ~exist('ignore_bounds', 'var')
        ignore_bounds = false;
    end
    [cBVI_to,...
     cBVI_opt_to,...
     cBVI_to_polypoints,...
     cBVI_opt_to_polypoints] = deal(cell(size(p.phi_def)));
    [to, to_opt] = calc_tlb_thirdorder(s,p,is_square);
    for shch = 1:length(p.cBVI)
        % skip if no cBVI (no third-order)
        if isempty(p.cBVI{shch}) || isempty(p.cBVI_opt{shch})
            continue;
        end
        % compute third-order effects (incl. exponential)
        cBVI_to{shch} = exp_map(p.cBVI{shch} + to{shch}{1});
        cBVI_opt_to{shch} = exp_map(p.cBVI_opt{shch} + to_opt{shch}{1});
        % error polygon
        cBVI_to_polypoints{shch} = repmat(p.cBVI{shch}, 1, 4) + to{shch}{2};
        cBVI_opt_to_polypoints{shch} = repmat(p.cBVI_opt{shch}, 1, 4) +...
                                       to_opt{shch}{2};
        % exponentiate polypoints
        for i = 1:4
            cBVI_to_polypoints{shch}(:,i) = exp_map(cBVI_to_polypoints{shch}(:,i));
            cBVI_opt_to_polypoints{shch}(:,i) = exp_map(cBVI_opt_to_polypoints{shch}(:,i));
        end

        % plot third-order effects
        hold on
        color = [234 14 30]/255;
        %color = [0, 183, 255]/255 + shch/length(p.cBVI)*([98, 0, 255] - [0, 183, 255])/255;
        if ~is_opt
            add_point_and_bounds(cBVI_to, cBVI_to_polypoints, shch, color, ~ignore_bounds);
        else
            add_point_and_bounds(cBVI_opt_to, cBVI_opt_to_polypoints, shch, color, ~ignore_bounds);
        end
        hold off
    end
end

function add_point_and_bounds(cBVI_to, cBVI_to_polypoints, shch, color, do_bounds)
    name = 'cBVI with Third Order ';
    if shch ~= 1
        name = [name num2str(shch)];
    end
    % TODO: remove the negative below if you have problems
    plot(cBVI_to{shch}(1), -cBVI_to{shch}(2), '+',...
         'MarkerSize', 15,...
         'LineWidth', 3,...
         'Color', color,...
         'DisplayName', name);
    if do_bounds
        try
            tmp = plot(polyshape(cBVI_to_polypoints{shch}(1,:),...
                                 cBVI_to_polypoints{shch}(2,:)),...
                        '-.', 'Color', color, 'HandleVisibility','off');
            alpha(tmp, 0.5);
        catch
            warning('Caught error due to 1D polyshape; plotting linear bounds');
            plot(cBVI_to_polypoints{shch}(1,[1 end]),...
                 cBVI_to_polypoints{shch}(2,[1 end]),...
                 '-_', 'Color', color, 'HandleVisibility','off');
        end
    end
end

function res = exp_map(vec)
    mat = [0 -vec(3) vec(1); vec(3) 0 vec(2); 0 0 0];
    res = mat_to_vec_SE2(expm(mat));
end