% creates output objects for dynamic visualization of error chains
function [joint_configurations,...
          expected_terms,...
          sep_error_terms,...
          error_points,...
          error_orders] = dynamic_error_chain(geometry,...
                                              configuration,...
                                              order,...
                                              truncate_to,...
                                              ranges)
    % inputs:
        % geometry:
            % elements: cell array of GroupElements encoding transforms
            % types: cell array of characters encoding transform types
                % 'j': joint
                % 'e': error
                % 'l': link
        % configuration: all scalar/sym arrays, multiplying their
        % respective transform types
            % links
            % joints
            % errors
        % order: integer order for both Taylor Series and BCH approx.
        % truncate_to: integer degree for term truncation
        % ranges: cell array for each symvar (in order), each element
        % containing the range for output plots
    % outputs:
        % joint_configurations: 
            % Nxdim array containing joint configurations
        % expected_terms:
            % cell array of polynomial terms in any supplied error 
            % variables, up to degree supplied by truncate_to
        % sep_error_terms:
            % cell array of terms equivalent to end_error_lge, but 
            % separated into quantities manifested as expected_terms
        % error_points: 
            % cell array, each entry containing plottable points for each 
            % error term
        % error_orders:
            % integer order of error term
        
    % generate points, error terms using error_chain.m
    [~,...
     ~,...
     expected_terms,...
     sep_error_terms,...
     joint_configurations] = error_chain(geometry,...
                                         configuration,...
                                         order,...
                                         truncate_to);
    % identify symvars
    err_vars = [];
    for t = 1:length(expected_terms)
        curr_vars = symvar(expected_terms{t});
        exist_vars = ismember(curr_vars, err_vars);
        err_vars = [err_vars curr_vars(~exist_vars)];
    end
    % confirm that ranges has supplied all error vars
    if length(err_vars) ~= length(ranges)
        error('dynamic_error_chain: variable ranges must be supplied for each symvar')
    end
    % preallocate error_points
    error_points = cell(size(expected_terms));
    error_orders = cell(size(expected_terms));
    % generate error points iteratively
    for t = 1:length(expected_terms)
        % pre-compute relevant values
        error_orders{t} = polynomialDegree(expected_terms{t});
        vars = symvar(expected_terms{t});
        var_idxs = find(ismember(vars, err_vars) == 1);
        num_vars = length(vars);
        % evaluate magnitude of local error term (using subs)
        vars_cell = num2cell(vars);
        range_cell = cell(size(vars_cell));
        for c = 1:length(vars)
            range_cell{c} = ranges{var_idxs(c)};
        end
        % range x dim matrix containing vectors
        term_vec = subs(sep_error_terms{t}, vars_cell, range_cell)';
        % rotate ea. vector into correct frame
        term_vec_rot = zeros(size(term_vec));
        for p = 1:size(term_vec, 1)
            point_ge = GroupElement(term_vec(p,:));
            ee_ge = GroupElement(joint_configurations(end,:));
            point_ge_rot = ee_ge.right_action(point_ge);
            [~, term_vec_rot(p,:)] = point_ge_rot.vector_true();
        end
        % assign error points
        error_points{t} = term_vec_rot;
    end
end