function vector = mat_to_vec_SE2_lie(matrix)
%Convert SE(2) matrices to a set of row vectors

vector = [matrix(1,3);matrix(2,3);matrix(2,1)];

end