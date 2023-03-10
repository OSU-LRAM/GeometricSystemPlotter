% does chain construction and error visualization using GroupElements
function [end_error_lae,...
          end_error_lge,...
          expected_terms,...
          sep_error_terms,...
          joint_configurations,...
          J] = error_chain(geometry, configuration, order, truncate_to)
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
    % outputs:
        % end_error_lae:
            % AlgebraElement encoding the joint error propogated to the end
            % effector
        % end_error_lge:
            % GroupElement after exponentiation of end_error_lae using
            % taylor series approximation to supplied order
        % expected_terms:
            % cell array of polynomial terms in any supplied error 
            % variables, up to degree supplied by truncate_to
        % sep_error_terms:
            % cell array of terms equivalent to end_error_lge, but 
            % separated into quantities manifested as expected_terms
        % joint_configurations:
            % matrix of nominal world space configurations of supplied
            % joints
        % J:
            % 3D array (pages of matrices) containing Jacobians of supplied
            % joints
    
    %% input sanitation
    % geometry checking
    if ~exist('geometry', 'var')
        error('error_chain: geometry must be provided');
    end
    if ~isfield(geometry, 'elements') || ~isfield(geometry, 'types')
        error('error_chain: geometry requires both elements and types');
    end
    if ~all(cellfun(@(A) isa(A, 'GroupElement'), geometry.elements))
        error('error_chain: geometry elements must be GroupElements')
    end
    if ~all(cellfun(@(A) isa(A, 'char'), geometry.types))
        error('error_chain: geometry types must be characters')
    end
    if ~exist('order', 'var')
        warning('error_chain: order not provided; assuming order 1');
        order = 1;
    end
    if ~exist('truncate_to', 'var')
        warning('error_chain: truncate_to not provided; truncated to degree 1');
        truncate_to = 1;
    end
    % create logical arrays for easy indexing
    link_logical = cellfun(@(c) strcmp(c,'l'), geometry.types);
    joint_logical = cellfun(@(c) strcmp(c,'j'), geometry.types);
    error_logical = cellfun(@(c) strcmp(c,'e'), geometry.types);
    % confirm existence of configuration
    if ~exist('configuration', 'var')
        configuration = struct;
    end
    if any(link_logical) && ~isfield(configuration, 'links')
        % assume links are as-given
        configuration.links = ones(1, sum(link_logical));
    end
    if any(joint_logical) && ~isfield(configuration, 'joints')
        % assume reference for joints
        configuration.joints = zeros(1, sum(joint_logical));
    end
    if any(error_logical) && ~isfield(configuration, 'errors')
        % assume reference for errors
        configuration.errors = zeros(1, sum(error_logical));
    end
    
    %% construction of Lie Algebra elements
    % get corresponding algebra elements from input group elements
    algebras = cell(1,length(geometry.elements));
    for i = 1:length(geometry.elements)
        % may be able to write as cellfun
        algebras{i} = geometry.elements{i}.to_AlgebraElement();
    end
    % apply scaling after log for cleanliness (is this admissible?)
    % this isn't very clean, but it's clear
    links = algebras(link_logical);
    joints = algebras(joint_logical);
    errors = algebras(error_logical);
    for l = 1:length(links)
        links{l} = AlgebraElement(configuration.links(l)*links{l}.vector);
    end
    for j = 1:length(joints)
        joints{j} = AlgebraElement(configuration.joints(j)*joints{j}.vector);
    end
    for e = 1:length(errors)
        errors{e} = AlgebraElement(configuration.errors(e)*errors{e}.vector);
    end
    algebras(link_logical) = links;
    algebras(joint_logical) = joints;
    algebras(error_logical) = errors;
    
    %% error math
    % check for error terms before jumping in
    if any(error_logical)
        %% symbolic math for error approximation
        % create symbolic error term
        dim = length(algebras{1}.vector);
        sigma_vec = sym('sigma', [1, dim]);
        assume(sigma_vec, 'real');
        sigma = AlgebraElement(sigma_vec);
        % use BCH to write LHS, RHS
        lhs_elements = [algebras(~error_logical) {sigma}];
        LHS = BCH_n_elements(lhs_elements, order);
        RHS = BCH_n_elements(algebras, order);
        % solve for error, sigma
        sigma_soln_struct = solve(LHS.vector == RHS.vector, sigma_vec);
        sigma_soln_vec = subs(sigma_vec, sigma_soln_struct);
        end_error_lae = AlgebraElement(sigma_soln_vec);
        % exponentiate
        end_error_lge = end_error_lae.to_GroupElement_Taylor(order);
        end_error_lge = end_error_lge.vector_small_angles();
        %% create polynomial of error terms (to search)
        % get supplied error variables vector
        err_vars = [];
        for e = 1:length(errors)
            curr_vars = symvar(errors{e}.matrix);
            exist_vars = ismember(curr_vars, err_vars);
            err_vars = [err_vars curr_vars(~exist_vars)];
        end
        % create possible polynomial terms
        full_poly = 1;
        for v = 1:length(err_vars)
            tmp_poly = 1;
            for d = 1:truncate_to
                tmp_poly = tmp_poly + err_vars(v)^d;
            end
            full_poly = full_poly*tmp_poly;
        end
        full_poly = expand(full_poly);
        % represent as ind. terms; truncate to desired degree
        expected_terms = children(full_poly);
        t = 1;
        while t <= length(expected_terms)
            if polynomialDegree(expected_terms{t}) > truncate_to
                expected_terms(t) = [];
            else
                t = t+1;
            end
        end
        %% search for expected polynomial terms
        % get terms in aggregate error matching expected terms
        end_error_poly = expand(end_error_lge.vector);
        end_error_terms = children(end_error_poly);
        sep_error_terms = cell(size(expected_terms));
        % searching for expected terms
        for term = 1:length(expected_terms)
            contribution = sym(zeros(size(end_error_terms,1),1));
            % search along ea. dimension
            for dim = 1:length(end_error_terms)
                % inspect ea. term of full end error polynomial
                for t = 1:length(end_error_terms{dim})
                    unit_quantity = end_error_terms{dim}{t}/expected_terms{term};
                    symvars = symvar(unit_quantity);
                    % if this is correct term, should be no err. terms after
                    % division
                    curr_exist_err_vars = ismember(symvars, err_vars);
                    if ~any(curr_exist_err_vars)
                        contribution(dim) = contribution(dim) + end_error_terms{dim}{t};
                    end
                end
            end
            sep_error_terms{term} = contribution;
        end
    else
        % no error terms; set outputs
        end_error_lae = AlgebraElement(zeros(size(joints{1}.vector)));
        end_error_lge = GroupElement(zeros(size(end_error_lae.vector)));
        expected_terms = 0;
        sep_error_terms = 0;
    end
    %% produce standard N_link_chain outputs (link configs, jacobians)
    % pre-allocate storage for joint configs, Jacobians
    joint_configurations = zeros(length(joints), length(joints{1}.vector));
    J = zeros(size(joint_configurations, 2),...
              size(joint_configurations,1),...
              size(joint_configurations,1));
    % keep track of current configuration (start as identity)
    current_config = GroupElement(zeros(1,size(joint_configurations, 2)));
    % count link/joint (for scalar indexing)
    j = 1;
    l = 1;
    % iterate thru group elements
    for e = 1:length(geometry.elements)
        % if error, skip
        if geometry.types{e} == 'e'
            continue
        % if link or joint, scale by supplied scalars
        elseif geometry.types{e} == 'l'
            trans_vec = geometry.elements{e}.vector*configuration.links(l);
            l = l + 1;
        elseif geometry.types{e} == 'j'
            trans_vec = geometry.elements{e}.vector*configuration.joints(j);
            % not incrementing j; will do later
        end
        % apply scaled transform to current config.
        transform = GroupElement(trans_vec);
        current_config = current_config.right_action(transform);
        % for joints, save relevant data
        if geometry.types{e} == 'j'
            % save configuration
            [current_config, conf_vec] = current_config.vector_true();
            joint_configurations(j,:) = conf_vec';
            % construct columns of J w.r.t base link
            % done vector-wise (not group-wise) because I have a slow brain
            for j_prev = 1:(j-1)
                diff_vec = conf_vec' - joint_configurations(j_prev,:);
                J(:, j_prev, j) = cross(diff_vec, joints{j_prev}.vector);
            end
            % increment joint index
            j = j + 1;
        end
    end
end