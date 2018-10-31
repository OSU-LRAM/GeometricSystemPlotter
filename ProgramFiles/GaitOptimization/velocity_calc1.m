function [L_i,L_is1,L_is2,L_ip1, dL_i_x,dL_is1_x,dL_i_y,dL_is1_y] = velocity_calc1(i,r,M,dM,shv,cons)
% This code create the local normalization


% Find the angle of horizental line
SE2 = @(t) [cos(t) -sin(t);sin(t) cos(t)];


% Fit the plane over the 3 points
xnew = [(r(i+1,:) - r(i-1,:))/(norm(r(i+1,:)-r(i-1,:))) 0];
v_needed = [(r(i+1,:) - r(i,:))/(norm(r(i+1,:) - r(i,:))) 0];

znew_temp = cross(xnew,v_needed);
znew = znew_temp/norm(znew_temp);

ynew_temp = cross(znew,xnew);
ynew = ynew_temp/norm(ynew_temp);

x0 = [1 0 0];
y0 = [0 1 0];
z0 = [0 0 1];
% Make the transformation matrix
Tr = [dot(x0,xnew) dot(x0,ynew) dot(x0,znew);
    dot(y0,xnew) dot(y0,ynew) dot(y0,znew);
    dot(z0,xnew) dot(z0,ynew) dot(z0,znew)];

% % Determine the angle from point 1 to point 3 for transforming the three
th = atan2(r(i+1,2)-r(i-1,2), r(i+1,1)-r(i-1,1));
Tr = SE2(th);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the length of Triangle sides
L_i = sqrt((r(i+1,:)-r(i,:))*((M{i}+M{i+1})/2)*(r(i+1,:)-r(i,:))');
L_is1 = sqrt((r(i,:)-r(i-1,:))*((M{i-1}+M{i})/2)*(r(i,:)-r(i-1,:))');
if i ~= 2
    L_is2 = sqrt((r(i-1,:)-r(i-2,:))*((M{i-1}+M{i-2})/2)*(r(i-1,:)-r(i-2,:))');
else
    L_is2 = 0;
end

    L_ip1 = sqrt((r(i+2,:)-r(i+1,:))*((M{i+2}+M{i+1})/2)*(r(i+2,:)-r(i+1,:))');


for j = 1:shv
    
    dr = zeros(1,shv);
    dr(j) = 1;
    
    dL_i{j} = 1*(-dr*((M{i}+M{i+1} )/2)*(r(i+1,:)-r(i,:))' + (r(i+1,:)-r(i,:))*((M{i}+M{i+1} )/2)*-dr' + (r(i+1,:)-r(i,:))*((dM{i,j} )/2)*(r(i+1,:)-r(i,:))')/(2*L_i);
    dL_is1{j} = 1*(dr*((M{i-1}+M{i} )/2)*(r(i,:)-r(i-1,:))' + (r(i,:)-r(i-1,:))*((M{i-1}+M{i} )/2)*dr' + (r(i,:)-r(i-1,:))*(( dM{i,j})/2)*(r(i,:)-r(i-1,:))')/(2*L_is1);

end
    

dL_i_x = cell2mat(dL_i)*Tr;
dL_is1_x = cell2mat(dL_is1)*Tr;
dL_i_x(2) = 0;
dL_is1_x(2) = 0;
dL_i_x = dL_i_x*inv(Tr);     % Velocity in tangent direction is zero
dL_is1_x = dL_is1_x*inv(Tr);     % Velocity only in tangent direction exist

dL_i_y = cell2mat(dL_i)*Tr;
dL_is1_y = cell2mat(dL_is1)*Tr;
dL_i_y(1) = 0;
dL_is1_y(1) = 0;
dL_i_y = dL_i_y*inv(Tr);     % Velocity in tangent direction is zero
dL_is1_y = dL_is1_y*inv(Tr);     % Velocity only in tangent direction exist


% if ~isempty(cons)
%     % First need to check the node
%     for j = 1:length(cons.path_points)
%         
%         if ismember(i,cons.path_points{j})
% 
%             T_rail_normalized = cons.T_rail{j}/norm(cons.T_rail{j});
%             dL1 = dot(dL1,T_rail_normalized)*T_rail_normalized;
%             dL2 = dot(dL2,T_rail_normalized)*T_rail_normalized;
%             
%             % Find the gradiant to keep the path on the rail
%             F = 10*cons.connecting_vector{j};
% 
%         else
% 
%             F = 0;
% 
%         end
%     
%     end
%     
% else
%     
%     F = 0;
%         
% end



end
