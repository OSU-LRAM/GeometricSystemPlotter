function matrix = vec_to_mat_SE2(vector)
%x-y-theta vector should be rows of the input matrix if multiple vector
%are supplied

%handle single column entry
if all(size(vector) == [3,1])
	vector = vector';
end

x = vector(:,1);
y = vector(:,2);
costheta = cos(vector(:,3));
sintheta = sin(vector(:,3));


%SE2 vector to matrix
matrix = zeros(3,3,size(x,1));

% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if isa(vector,'sym')
    matrix = sym(matrix);
end

matrix(1,1,:) = costheta;
matrix(1,2,:) = -sintheta;
matrix(1,3,:) = x;

matrix(2,1,:) = sintheta;
matrix(2,2,:) = costheta;
matrix(2,3,:) = y;

matrix(3,3,:) = ones(1,1,size(x,1));


end
