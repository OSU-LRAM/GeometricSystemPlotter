function matrix = vec_to_mat_SE2_lie(vector)

theta = vector(3);
x = vector(1);
y = vector(2);

matrix = [0,-theta,x;theta,0,y;0,0,0];

end
