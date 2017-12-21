% Empirical derivation of local connection matrix
% requires pre-existing body velocity and shape velocity variables

function output = empirical_local_connection(xd,a,ad,num_grid,min_max,point_range)

% xd ------- body velocity in the x direction. nx1 array
% a -------- shape variable inputs. nx2 array
% ad ------- shape velocity inputs. nx2 array
% num_grid - number of grid points. scalar value
% min_max -- minimum and maximum for grid size. 1x2 vector

% replace this with actual data later
% n = 100;
n = size(xd,1);
% m = 2;
% m = size(a,2);

% change which body velocity components are evaluated for x, y, th
% xd = zeros(n,1); % matrix of sampled body velocities, 1xn, x direction. replace with actual data.
% a = zeros(n,2);  % shape coordinate matrix, two dimenstions (a1 and a2)
% ad = zeros(n,2); % shape velocity matrix, two directions (a1d and a2d)

% generate deltas given grid and tolerance
ng = num_grid; % number of grid spaces
agrid = linspace(min_max(1),min_max(2),ng); % generate grid along one dimension

[a1,a2] = ndgrid(agrid,agrid);

srange = point_range;

% initialize deltas for later use
d1 = [];
d2 = [];

A = zeros(ng,ng);

for i = 1:ng
    
    for j = 1:ng
        
        pt = [a1(i,j) a2(i,j)]; % 1x2 coordinate matrix for grid location
        
        for k = 1:n % find all points within boundary of
            
            % if the data point is within the defined range.
            % uses circle as boundary; maybe switch out to grid? multiple
            % instances of points?
%             if ((pt(1) - a(1))^2 + (pt(2) - a(2))^2) <= srange
            if abs(pt(1) - x(k)) <= srange && abs(pt(2) - y(k)) <= srange
                
                % add [x_index,delta] to end matrix
                d1 = [d1;k,(a(1)-pt(1))]; %#ok<AGROW> 
                d2 = [d2;k,(a(2)-pt(2))]; %#ok<AGROW> 
                
            end
        end
        ad_matrix = [];
        for k2 = 1:size(d1,1)
            
            idx = d1(k2); % find index which xd array uses
            xbd(k2,1) = xd(idx); % find associated body velocity value
            ad_matrix = [ad_matrix;... % build ad matrix for pseudo-inverse
                ad(idx,1),ad(idx,2),ad(idx,1)*d1(k2,2),ad(idx,1)*d2(k2,2),ad(idx,2)*d1(k2,2),ad(idx,2)*d2(k2,2)];
            %   ad1       ad2       ad1*d1            ad1*d2            ad2*d1            ad2*d2
            
        end
        
        A_mat = pinv(ad_matrix)*xbd;
        A1ad2 = A_mat(4); % select dA1/da2 and dA2/da1
        A2ad1 = A_mat(5);
        
        % here is where the important values are changed. maybe use more
        % data in future?
        A(i,j) = A2ad1 - A1ad2; % surface value for curl + partial
        
        
        
    end
end

% output A as a num_grid x num_grid array to be used as surface later. 
output = A;

end



