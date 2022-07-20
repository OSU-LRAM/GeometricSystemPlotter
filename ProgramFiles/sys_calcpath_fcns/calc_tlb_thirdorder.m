% calculates approximation of third-order gait effects
% returns third order effects in cell:
    % third_order{shch} refers to effects for given shape change
    % third_order{shch}{1} provides nominal third-order effect vector
    % third_order{shch}{2} provides all worst-case third-order effect
    %                      vectors (shaped into a matrix for plotting)
function [third_order, third_order_opt] = calc_tlb_thirdorder(s,p,is_square)
    % sanitize
    if ~exist('is_square', 'var')
        is_square = false;
    end
    % setup
    [third_order, third_order_opt] = deal(cell(size(p.phi_def)));
    % iterate thru shape changes
    for shch = 1:length(third_order)
        % skip if no cBVI
        if isempty(p.cBVI{shch}) || isempty(p.cBVI_opt{shch})
            warning('No cBVI available for gait, skipping third order calculation');
            continue;
        end
        % get unit tangent from beginning of gait (FD approx.)
        first_two = p.phi_fun_full{shch}(p.time_full{shch}(1:2))';
        init_vel = (first_two(:,2) - first_two(:,1));
        % initial step * unit step velocity * time for 1st quarter
        steps = length(p.time_full{shch});
        init_vel = init_vel * steps / p.time_full{shch}(end) * p.time_full{shch}(floor(steps/4));

        % find center of gait
        alphas = p.phi_fun_full{shch}(p.time_full{shch});
        center = mean(alphas, 1);
        
        % local connection estimates
        LC_center = {zeros(3,2), zeros(3,2)};
        A = {s.vecfield.eval.content.Avec,...;
             s.vecfield.eval.content.Avec_optimized};
        % evaluate LC at center (mean) of gait, in orig. and opt. coords
        % interp assumes two shape vars
        for shvar = 1:length(init_vel)
            for dim = 1:3
                for i = 1:2 %orig, opt. coords
                    LC_interp = interp2(s.grid.eval{2}, s.grid.eval{1},...
                                        A{i}{dim, shvar},...
                                        center(1), center(2));
                    LC_center{i}(dim, shvar) = LC_center{i}(dim, shvar) +...
                                               LC_interp;
                end
            end
        end
        
        % use mean of entire LC, communicating magnitude and direction
        % jury's out on whether using the entire LC is a good idea
        %for i = 1:numel(A{1})
        %    for c = 1:2 %orig, opt. coords
        %        LC_center{c}(i) = mean(A{c}{i}, 'all');
        %    end
        %end

        % tools for approx.
        R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
        alpha = cell(1,2);
        beta = alpha;
        % do third order approx.
        for i = 1:2 %orig, opt. coords
            if is_square
                alpha{i} = LC_center{i} * init_vel; %init. trajectory
                beta{i} = LC_center{i} * R(pi/2) * init_vel; %second leg
            else
                % evaluate LC at unit TV rotated by 1/8 and 3/8 (approx. alpha, beta)
                % scaled by pi/4, mapping square gait to circle of equal diameter
                alpha{i} = pi/4 * LC_center{i} * R(pi/4) * init_vel;
                beta{i} = pi/4 * LC_center{i} * R(3*pi/4) * init_vel;
            end
        end

        % lie bracket definition
        lb = @(x,y) [y(3)*x(2) - x(3)*y(2); x(3)*y(1) - y(3)*x(1); 0];
        % generate storage
        [third_order{shch}, third_order_opt{shch}] = deal(cell(1,2));
        % nominal location
        third_order{shch}{1} = lb(alpha{1}+beta{1}, p.cBVI{shch})/2;
        third_order_opt{shch}{1} = lb(alpha{2}+beta{2}, p.cBVI_opt{shch})/2;
        % worst-case
        % lie bracket definition (worst-case)
        lb = @(x,y) [abs(y(3)*x(2)) + abs(x(3)*y(2)); abs(x(3)*y(1)) + abs(y(3)*x(1)); 0];
        pos_worst = lb(abs(alpha{1}) + abs(beta{1}), p.cBVI{shch})/2;
        pos_worst_opt = lb(abs(alpha{2}) + abs(beta{2}), p.cBVI_opt{shch})/2;
        % combinations of worst-case in +/- x,y only
        combinations = [eye(3); ...
                        [-1 0 0; 0 1 0; 0 0 1];...
                        [1 0 0; 0 -1 0; 0 0 1];...
                        -eye(3)];
        third_order{shch}{2} = reshape(combinations*pos_worst, [3,4]);
        third_order_opt{shch}{2} = reshape(combinations*pos_worst_opt, [3,4]);
    end
end