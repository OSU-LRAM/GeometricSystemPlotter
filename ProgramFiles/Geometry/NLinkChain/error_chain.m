% does chain construction and error visualization using GroupElements
function end_error_lae = error_chain(geometry, configuration, order)
    % inputs:
        % geometry:
            % elements: cell array of GroupElements encoding transforms
            % types: cell array of characters encoding transform types
                % 'j': joint
                % 'e': error
                % 'l': link (can really be anything)
        % configuration: all scalar/sym arrays, multiplying their
        % respective transform types
            % links
            % joints
            % errors
        % order: integer order for both Taylor Series and BCH approx.
    
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
    % create logical arrays for easy indexing
    link_logical = cellfun(@(c) strcmp(c,'l'), geometry.types);
    joint_logical = cellfun(@(c) strcmp(c,'j'), geometry.types);
    error_logical = cellfun(@(c) strcmp(c,'e'), geometry.types);
    % confirm existence of configuration
    if ~exist('configuration', 'var')
        configuration = struct;
    end
    if any(link_logical) && ~isfield(configuration, 'links')
        configuration.links = ones(1, sum(link_logical));
    end
    if any(joint_logical) && ~isfield(configuration, 'joints')
        configuration.joints = ones(1, sum(joint_logical));
    end
    if any(error_logical) && ~isfield(configuration, 'errors')
        configuration.errors = ones(1, sum(error_logical));
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
    
    %% symbolic math for error approximation
    % create symbolic error term
    dim = length(algebras{1}.vector);
    sigma_vec = sym('sigma', [1, dim]);
    sigma = AlgebraElement(sigma_vec);
    % use BCH to write LHS, RHS
    lhs_elements = [algebras(~error_logical) {sigma}];
    LHS = BCH_n_elements(lhs_elements, order);
    RHS = BCH_n_elements(algebras, order);
    % solve for error, sigma
    sigma_soln_struct = solve(LHS.vector == RHS.vector, sigma_vec);
    sigma_soln_vec = subs(sigma_vec, sigma_soln_struct);
    end_error_lae = AlgebraElement(sigma_soln_vec);
    % exponentiate with T.S. approx.
    end_error_lge = end_error_lae.to_GroupElement_Taylor(order);
    
    %% break solution apart in terms of errors
    % break into terms: linear, bilinear, quadratic... (in ea. err)
    % plot contributions (in an interactable plot environment?)
    
    % existing outputs:
        % link configurations
        % Jacobians of each link
    % new outputs:
        % error term(s)?
        % plotted errors
end