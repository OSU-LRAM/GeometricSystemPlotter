function dMdq = partial_mass_matrix(J_full,dJdq,local_inertias)
% Calculates the partial of the full mass matrix M with respect to the
% configuration variables q = [g_circ r]', where g_circ is the system's
% global coordinate configuration and r is the local shape-space
% configuration.
%
% Inputs for a system with m links and q of size (n by 1):
%
%   J_full: Full Jacobian matrix for the system, including entries for all
%       configuration variables q; (m by 1) cell array where each cell contains
%       a (3 by n) matrix representing the Jacobian of link m.
%
%   dJdq: Partial of full Jacobian with respect to the configuration
%       variables q; (m by 1) cell array where each entry corresponds to
%       the partial derivative of the Jacobian for link m, which is another
%       cell array of size (n by n) where each entry (j,k) corresponds to a
%       three-element vector that represents the partial of Jacobian column
%       j with respect to configuration variable k.
%
%   local_inertias: The inertia tensor of the link as measured in its fixed
%       coordinate frame, which includes the added mass from the surrounding
%       fluid. Cell array where the ith cell corresponds to the inertia
%       matrix for link i.
%
% Output for a system with q of size (n by 1):
%
%   dMdq: Partial of full mass matrix with respect to the configuration
%       variables q; (n by 1) cell array where the ith cell contains an (n by n)
%       matrix corresponding to dMdq_i, or the partial of the mass matrix with
%       respect to the ith configuration variable in q

% Number of links is equal to the number of Jacobian matrices in J_full, as
% it contains one matrix for each link
num_links = length(J_full);
% Number of configuration variables q
num_q = size(dJdq{1},1);

% We need a partial mass matrix for each configuration variable in q, so
% create a temporary prototype with the appropriately sized matrix to enter
% into preallocated space
dMtemp = zeros(num_q,num_q);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or(isa(J_full{1},'sym'),isa(local_inertias{1},'sym'))
    dMtemp = sym(dMtemp);
    cell2func = @cell2sym;
else
    cell2func = @cell2mat;
end
% Now fill the preallocated memory with the temporary prototypes
dMdq = cell(num_q,1); dMdq(1:end) = {dMtemp};

% The mass matrix for a mobile system is given by M = J'*mu*J, where J is
% the full Jacobian and mu is the local inertia tensor, summed across each 
% of the links. We can obtain dMdq simply by performing the chain rule to 
% obtain dMdq = dJdq'*mu*J + J'*mu*dJdq, summed across all of the links,
% and performed for each configuration variable q, as done below.
for q = 1:num_q
    for link = 1:num_links
        dJtemp = cell2func(dJdq{link}(:,q)');
         dMdq{q} = dMdq{q} + dJtemp' * local_inertias{link} * J_full{link} ...
            + J_full{link}' * local_inertias{link} * dJtemp;
    end
    if isa(dMtemp,'sym')
        dMdq{q} = simplify(dMdq{q},'Steps',10);
    end
end