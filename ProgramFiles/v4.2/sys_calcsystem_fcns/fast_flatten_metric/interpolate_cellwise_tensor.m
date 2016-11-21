function T = interpolate_cellwise_tensor(x,y,x_p,y_p,T)
% Take a tensor defined as a cell array with the values in each index, and
% return its value at a given point.

T = cell2mat(cellfun(@(T_i) interpn(x,y,T_i,x_p,y_p),T,'UniformOutput',false));

end