% generate shape functions for an N-dimensional hypercubic element

function [shape_functions, shape_dfunctions] = hypercube_element_shape_functions(N)

	% prime arrays of shape functions and dfunctions
	shape_functions = cell(2^N,1);
	shape_dfunctions = cell(2^N,N);
	
	% set up node orderings for the shape functions
    % uses similar binary setup as in hypercube_mesh.m
    bin_strs = dec2bin(0:2^N - 1);
    bin_mat = zeros(size(bin_strs));
    bin_mat(:) = arrayfun(@str2num, bin_strs(:));
    ordering = fliplr(2*bin_mat-1);
	
	% convert ordering into a cell array for processing inside the shape
	% functions
	ordering = num2cell(ordering);
	
	% construct the shape functions and dfunctions
	for i = 1:length(shape_functions)
		
		shape_functions{i} = @(X) 1/(2^N) * shape_fun_prim(ordering(i,1:N),X);
		
		for j = 1:size(shape_dfunctions,2)
			
			shape_dfunctions{i,j} = @(X) 1/(2^N) * shape_dfun_prim(ordering(i,1:N),X,j);
		
		end
	end

end


function psi = shape_fun_prim(signs,X)
% (n-)linear shape function for a hypercubic element with the input signs
% giving the signs for the polynomial
	
	% Make sure that there are exactly as many inputs as signs
	if numel(signs) ~= numel(X)
		error(['This shape fun was generated over a ' num2str(numel(signs)) ' space, but has been given a ' num2str(numel(X)) '-dimensional input'])
	end

	% Put the inputs and signs into the shape function
	psi_factor = cellfun(@(x,y) 1 + (x*y),(reshape(signs,size(X))),X,'UniformOutput',false);
	
	%%%%%%%%
	% Calculate the value of the shape function at each input
	
	% Concatenate the inputs along a higher dimension than they contain
	% already
	psi_stack = cat(numel(size(psi_factor{1}))+1,psi_factor{:});
	
	% Multiply down
	if numel(X)~=1
		psi = prod(psi_stack,numel(size(psi_stack)));
	else
		psi = psi_stack;
	end
	
end

function dpsi = shape_dfun_prim(signs,X,dcoord)
% Derivative of an (n-)linear shape function for a hypercubic element, with
% respect to the coordinate enumerated as dcoord
	
	% Make sure that there are exactly as many inputs as signs
	if numel(signs) ~= numel(X)
		error(['This shape derivative fun was generated over a ' num2str(numel(signs)) ' space, but has been given a ' num2str(numel(X)) '-dimensional input'])
	end
		
	% Put the inputs and signs into the shape derivative function
	dpsi_factor = cellfun(@(x,y) 1 + (x*y),(reshape(signs,size(X))),X,'UniformOutput',false);
	
	% Replace the factor corresponding to the dcoord with just the sign for
	% that term (the effect of taking a derivative with respect to that
	% shape function, given the linear-polynomial definition of the
	% function
	dpsi_factor{dcoord} = signs{dcoord}*ones(size(X{dcoord}));
	
	%%%%%%%%
	% Calculate the value of the shape derivative function at each input
	
	% Concatenate the inputs along a higher dimension than they contain
	% already
	dpsi_stack = cat(numel(size(dpsi_factor{1}))+1,dpsi_factor{:});
	
	% Multiply down
	if numel(X)~=1
		dpsi = prod(dpsi_stack,numel(size(dpsi_stack)));
	else
		dpsi = dpsi_stack;
	end
	
end