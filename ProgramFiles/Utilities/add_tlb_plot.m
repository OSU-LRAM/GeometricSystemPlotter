function [cBVI_to, cBVI_opt_to] = add_tlb_plot(s,p,is_opt,is_square)
    % sanitation
    if ~ishandle(1)
        error('Open trajectory graph before adding third-order effects');
    elseif ~exist('is_opt','var')
        is_opt = false;
    elseif ~exist('is_square', 'var')
        is_square = false;
    end
    [cBVI_to,...
     cBVI_opt_to,...
     cBVI_to_polypoints,...
     cBVI_opt_to_polypoints] = deal(cell(size(p.phi_def)));
    for shch = 1:length(p.cBVI)
        % compute third-order effects
        [to, to_opt] = calc_tlb_thirdorder(s,p,is_square);
        cBVI_to{shch} = p.cBVI{shch} + to{shch}{1};
        cBVI_opt_to{shch} = p.cBVI_opt{shch} + to_opt{shch}{1};
        cBVI_to_polypoints{shch} = repmat(p.cBVI{shch}, 1, 4) + to{shch}{2};
        cBVI_opt_to_polypoints{shch} = repmat(p.cBVI_opt{shch}, 1, 4) +...
                                       to_opt{shch}{2};

        % plot third-order effects
        hold on
        color = [0, 183, 255]/255 + shch/length(p.cBVI)*([98, 0, 255] - [0, 183, 255])/255;
        if ~is_opt
            plot(cBVI_to{shch}(1), cBVI_to{shch}(2), 'o',...
                'MarkerSize', 7,...
                'LineWidth', 2,...
                'Color', color,...
                'DisplayName', ['cBVI_{to}' num2str(shch)]);
            try
                tmp = plot(polyshape(cBVI_to_polypoints{shch}(1,:),...
                                     cBVI_to_polypoints{shch}(2,:)),...
                            '-.', 'Color', color);
                alpha(tmp, 0.5);
            catch
                warning('Caught error due to empty polyshape; skipping plot');
            end
        else
            plot(cBVI_opt_to{shch}(1), cBVI_opt_to{shch}(2), 'o',...
                'MarkerSize', 7,...
                'LineWidth', 2,...
                'Color', color,...
                'DisplayName', ['cBVI_{opt,to}' num2str(shch)]);
            try
                tmp = plot(polyshape(cBVI_opt_to_polypoints{shch}(1,:),...
                                     cBVI_opt_to_polypoints{shch}(2,:)),...
                            '-.', 'Color', color);
                alpha(tmp, 0.5);
            catch
                warning('Caught error due to empty polyshape; skipping plot');
            end
        end
        hold off
    end
end