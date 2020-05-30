function [A_out, h, J, J_full, omega] = Nonholonomic_connection_discrete(geometry,jointangles,numden)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields required for this function are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are
%           specified in the documentation of N_link_chain.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%        constraint_list: Structure array containing information about where the
%            constraints on the system are located and how they are oriented.
%            One structure array element for each constraint
%  
%            Fields in each array element are:
%  
%            link_number: link on which the constraint is located
%  
%            constraint_direction: direction of motion (in the link's body
%                coordinates) that is prohibited by the constraint
%
%   jointangles: Angles of the joints between the links, defining the
%       chain's current shape.
%
%   numden: Specification of whether the function should return the
%       numerator or denominator of the local connection, or if these should
%       be combined and the full expression for the local connection
%       returned.
%
% Outputs:
%
%   A_out : The "local connection" matrix, or "Locomotion Jacobian". This matrix
%       maps joint angular velocities ("shape velocities") to body
%       velocities of the system's base frame as 
%
%            g_b = -A * alphadot
%
%       with g_b the body velocity and alphadot the joint angular velocity.
%
%       Because there may be singularities in the connection, the input
%       argument 'numden' allows the user to specify that A be split into a
%       numerator and a denominator, with A = A_num ./ A_den, and only one
%       of these pieces returned; this allows the 'smart_divider' function
%       in sysplotter to identify and special-case the singularities.
%
%       (Note the negative sign in this equation; sysplotter code uses the
%       classical formulation of the geometric mechanics equations, which
%       includes this negative sign.
%
%   [h,J,J_full]: location and Jacobians for the links on the chain. These
%       are passed through from N_link_chain; see the documentation on that
%       function for more details
%
%   omega: The "Pfaffian constraint matrix" from which the local connection
%       A is calculated. This matrix is a linear map from system body and
%       shape velocities to net external forces acting on the base frame,
%       which must be zero for all achievable motions of the system


    % Extract the constraint list from the geometry
    constraint_list = geometry.constraint_list;

    % Default to providing the full connection without splitting it into
    % numerator and denominator
    if ~exist('numden','var')
        numden = 'combined';
    end

    %%%%
    % First, get the positions of the links in the chain and their
    % Jacobians with respect to the system parameters
    if isfield(geometry,'subchains')
        [h, J, J_full] = branched_chain(geometry,jointangles);
    else
        [h, J, J_full] = N_link_chain(geometry,jointangles);
    end

    
    %%%%
    % Next, iterate over the first three constraints. For each constraint,
    % multiply the local direction whose motion is prohibited by the
    % constraint into the full body-velocity Jacobian for the link; this
    % gives the contribution of the local (physical) constraint to the
    % constraints on the generalized coordinates
    
    % Preallocate Pfaffian constraint.
    % Pfaffian constraint will be of the same dimensionas the full link Jacobians
    omega = zeros(size(J_full{1}));
    
    % Prevent Matlab from playing tricks with imaginary numbers on symbolic
    % inputs and from complaining about assumptions on constants
    if or( isa(jointangles,'sym'), isa(geometry.linklengths,'sym') )
        omega = sym(omega);
        symjointangles = symvar(jointangles);
        if ~isempty(symjointangles)
            assume(symjointangles,'real')
        end
        symlinklengths = symvar(geometry.linklengths);
        if ~isempty(symlinklengths)
            assume(symlinklengths,'real')
        end

        warning('off','symbolic:sym:sym:AssumptionsOnConstantsIgnored')
    end
    
    % Iterate over the constraints, filling rows of the Pfaffian from
    % projections of the Jacobians onto the constraint
    for idx = 1:3
        omega(idx,:) = torow(constraint_list(idx).constraint_direction) * J_full{constraint_list(idx).link_number};
    end
    
    % Prevent symbolic expression from growing too cumbersome
    if isa(omega,'sym')
        omega = simplify(omega);
    end
    
    %%%%
    % Now separate out the the pfaffian into position and shape blocks
    omega_g = omega(:,1:3);
    omega_alpha = omega(:,4:end);
    
    %%%%
    % Local connection is A = omega_g \ omega_alpha. When singularities are
    % present in the constraints, they show up as omega_g not being
    % invertible, which corresponds to det(omega_g) going to zero. To
    % handle this numerically, we will calculate det(omega_g) separately
    % from the rest of inv(omega_g).
    
    % Calculate the determinant of the position block of omega
    det_omega_g = det(omega_g);
    
    %%%
    % Calculate the rest of inv(omega_g)
    
    % preallocate a matrix that is the same size as omega_g
    inv_omega_g_scaled = omega_g;
    
    % iterate over the entries of this inverse-like matrix. In each entry,
    % place the determinant of the submatrix of omega_g that doesn't
    % include the element we are filling
    for idx = 1:3
        for idx2 = 1:3
            inv_omega_g_scaled(idx,idx2) = subdet(omega_g,idx,idx2);
        end
    end
    
    
    % To get the numerator of the connection, multiply the scaled-inverse
    % of omega_g by omega_alpha
    A_num = inv_omega_g_scaled * omega_alpha;
    
    % To get the denominator of the connection, multiply det_omega_g by a
    % matrix of ones that is the same size as the numerator of the
    % connection
    A_den = det_omega_g * ones(size(A_num));

    % Simplify numerator and denominator
    if isa(omega,'sym')
        A_num = simplify(A_num,'Steps',50);
        A_den = simplify(A_den,'Steps',30);
    end
    
    % Return, the numerator, the denominator, or their combination
    switch numden
        case 'combined'
            A_out = A_den.\A_num;
        case 'num'
            A_out = A_num;
        case 'den'
            A_out = A_den;
    end

end

function subdetM = subdet(M,a,b)
% Find the determinant of the submatrix of 3x3 matrix M that does not include indices
% a and b

% Get the positively-ordered pair of indices that are not the input
% indices. Positively-ordering the pairs means that we don't need to monkey
% around with signs-of-cofactors
w = mod((b-1)+1,3)+1;
x = mod((b-1)+2,3)+1;

y = mod((a-1)+1,3)+1;
z = mod((a-1)+2,3)+1;

% Get the submatrix that corresponds to the indices selected. Note that the
% columns selected are the rows not present, and vice versa; this is
% because of the transpose operation in the inverse
subM = M([w x],[y z]);

% Get the determinant of this submatrix
subdetM = det(subM);

    % Prevent symbolic expression from growing too cumbersome
    if isa(subdetM,'sym')
        subdetM = simplify(subdetM);
        warning('off','symbolic:sym:sym:AssumptionsOnConstantsIgnored')

    end

end


