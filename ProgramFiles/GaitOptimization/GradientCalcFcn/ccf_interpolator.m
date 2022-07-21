function ccf = ccf_interpolator(y,s,n,dimension,direction)
% Preallocating memory for variables which we will need in further
% calculation 
yvalues=cell(n,dimension); % Cell representation of the coordinates of all points forming the gait
interpstateccf=cell(1,dimension); % Variable which will store the ccf function grid used for interpolation
ccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store ccf function at each point

% Interpolation to calculate all the variables needed for gradient
% calculation
for i=1:1:n
    for j=1:1:dimension
        yvalues{i,j}=y(i,j);
    end
end

y_for_interp = mat2cell(y,size(y,1),ones(1,size(y,2)));

for j=1:1:dimension
    interpstateccf{j}=s.grid.eval{j,1};
end

for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y_for_interp{:},'spline');
end


end