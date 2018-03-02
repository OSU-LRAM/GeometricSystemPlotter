function matrix = vec_to_mat_se2(vector)
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


%se2 vector to matrix
matrix = zeros(3,3,size(x,1));

matrix(1,1,:) = costheta;
matrix(1,2,:) = -sintheta;
matrix(1,3,:) = x;

matrix(2,1,:) = sintheta;
matrix(2,2,:) = costheta;
matrix(2,3,:) = y;

matrix(3,3,:) = ones(1,1,size(x,1));

%     matrix = [cos(vector(3)) -sin(vector(3)) vector(1);
%         sin(vector(3)) cos(vector(3)) vector(2);
%         0   0   1];
end
