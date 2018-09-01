function vector = mat_to_vec_SE2(matrix)
%Convert SE(2) matrices to a set of row vectors

x = permute(matrix(1,3,:),[3 2 1]);
y = permute(matrix(2,3,:),[3 2 1]);
theta = atan2(permute(matrix(2,1,:),[3 2 1]),permute(matrix(1,1,:),[3 2 1]));

vector = [x y theta];