function [kappa_final1,kappa_final2,F, dkappa1,dkappa2] = curvature_calc4(i,r,M,P,cons,Gamma,dM) %-----Removed method as fifth input

a = P.alpha_field{1};
b = P.alpha_field{2};

Mp_raw = P.M_t;

% This code create the local normalization

% Find the angle of horizental line
SE2z = @(t) [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];
SE2y = @(t) [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
SE2x = @(t) [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];

SE2 = @(t) [cos(t) -sin(t);sin(t) cos(t)];

% Singular Value Decomposition
[u,s,v] = svd(M{i});

% Tissot Indicatrix
T = u*s^-0.5*v';

% Transform the metrix with Tissot

% Transform the coordinate with Tissot
U_temp = r*inv(T);

U_tempp1 = r*inv(T);

U_temps1 = r*inv(T);

th = atan2(U_temp(i+1,2)-U_temp(i-1,2), U_temp(i+1,1)-U_temp(i-1,1));

thp1 = atan2(U_tempp1(i+2,2)-U_tempp1(i,2), U_tempp1(i+2,1)-U_tempp1(i,1));

ths1 = atan2(U_temps1(i,2)-U_temps1(i-2,2), U_temps1(i,1)-U_temps1(i-2,1));

th1 = atan2(r(i+1,2)-r(i-1,2), r(i+1,1)-r(i-1,1));
th1p1 = atan2(r(i+2,2)-r(i,2), r(i+2,1)-r(i,1));
th1s1 = atan2(r(i,2)-r(i-2,2), r(i,1)-r(i-2,1));

% Fit the plane over the 3 points
% xnew = [(U_temp(i+1,:) - U_temp(i-1,:))/(norm(U_temp(i+1,:)-U_temp(i-1,:))) 0];
% v_needed = [(U_temp(i+1,:) - U_temp(i,:))/(norm(U_temp(i+1,:) - U_temp(i,:))) 0];
% 
% znew_temp = cross(xnew,v_needed);
% znew = znew_temp/norm(znew_temp);
% 
% ynew_temp = cross(znew,xnew);
% ynew = ynew_temp/norm(ynew_temp);
% 
% if isnan(znew)
%     znew = zeros(3,1);
%     ynew = zeros(3,1);
% end
%     
% 
% x0 = [1 0 0];
% y0 = [0 1 0];
% z0 = [0 0 1];
% % Make the transformation matrix
% Tr = [dot(x0,xnew) dot(x0,ynew) dot(x0,znew);
%     dot(y0,xnew) dot(y0,ynew) dot(y0,znew);
%     dot(z0,xnew) dot(z0,ynew) dot(z0,znew)];
% 
% Tr = Tr(1:2,1:2);

Tr = SE2(th);

Tr1 = SE2(th1);
Tr1p1 = SE2(th1p1);
Tr1s1 = SE2(th1s1);

Trp1 = SE2(thp1);

Trs1 = SE2(ths1);

% Transform the coordinate with the above rotational angle
U = U_temp*Tr;

Up1 = U_tempp1*Trp1;

Us1 = U_temps1*Trs1;

U1 = r*(SE2(th1));
U1p1 = r*(SE2(th1p1));
U1s1 = r*(SE2(th1s1));

%-----Pretty sure Hossein said that we'll only use Numerical case-----%
% switch method
%     
%     case 'Analytical'   % If you have your metric equation you can use this option
%         
%         % Analytical function goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'Numerical'    % When the metric is given for all the configuration you can use this option

    Gamma2 = Gamma{i};
    Gammap1 = Gamma{i+1};
    Gammas1 = Gamma{i-1};

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain the length of Triangle sides
L_i = sqrt((r(i+1,:)-r(i,:))*((M{i}+M{i+1})/2)*(r(i+1,:)-r(i,:))');
L_is1 = sqrt((r(i,:)-r(i-1,:))*((M{i-1}+M{i})/2)*(r(i,:)-r(i-1,:))');

L_is2 = sqrt((r(i-1,:)-r(i-2,:))*((M{i-1}+M{i-2})/2)*(r(i-1,:)-r(i-2,:))');
L_ip1 = sqrt((r(i+2,:)-r(i+1,:))*((M{i+2}+M{i+1})/2)*(r(i+2,:)-r(i+1,:))');


% Speed squared required for finding Geodesic Curvature
v_it = (0.5*L_i + 0.5*L_is1)^2; 

% Find the change of angle
y1 = U(i,2) - U(i-1,2);
y2 = U(i+1,2) - U(i,2);
x1 = U(i,1) - U(i-1,1);
x2 = U(i+1,1) - U(i,1);


% next point
y1p1 = Up1(i+1,2) - Up1(i,2);
y2p1 = Up1(i+2,2) - Up1(i+1,2);
x1p1 = Up1(i+1,1) - Up1(i,1);
x2p1 = Up1(i+2,1) - Up1(i+1,1);


% % Previous point
y1s1 = Us1(i-1,2) - Us1(i-2,2);
y2s1 = Us1(i,2) - Us1(i-1,2);
x1s1 = Us1(i-1,1) - Us1(i-2,1);
x2s1 = Us1(i,1) - Us1(i-1,1);



%% Find the middle point of triangle in main coordinate (U1)
Mid_pts = [U1(i,1) (U1(i+1,2)+U1(i-1,2))/2];
Mid_ptsp1 = [U1p1(i+1,1) (U1p1(i+2,2)+U1p1(i,2))/2];
Mid_ptss1 = [U1s1(i-1,1) (U1s1(i,2)+U1s1(i-2,2))/2];

% Take the point to original coordinate to find the metric
Mid_pts_org = Mid_pts*inv(Tr1);
Mid_pts_orgp1 = Mid_ptsp1*inv(Tr1p1);
Mid_pts_orgs1 = Mid_ptss1*inv(Tr1s1);

% Find the metric
M_mid_org = Granular_metric_calc2({Mid_pts_org(1),Mid_pts_org(2)},Mp_raw,a,b);
M_mid_orgp1 = Granular_metric_calc2({Mid_pts_orgp1(1),Mid_pts_orgp1(2)},Mp_raw,a,b);
M_mid_orgs1 = Granular_metric_calc2({Mid_pts_orgs1(1),Mid_pts_orgs1(2)},Mp_raw,a,b);

% Transform the metric to rotated coordinate
M_mid = inv(SE2(th1))*M_mid_org*inv(SE2(th1))';
M1_i = inv(SE2(th1))*M{i}*inv(SE2(th1))';
M1_ip1 = inv(SE2(th1))*M{i+1}*inv(SE2(th1))';
M1_is1 = inv(SE2(th1))*M{i-1}*inv(SE2(th1))';

M_midp1 = inv(SE2(th1p1))*M_mid_orgp1*inv(SE2(th1p1))';
M1p1_i = inv(SE2(th1p1))*M{i}*inv(SE2(th1p1))';
M1p1_ip1 = inv(SE2(th1p1))*M{i+1}*inv(SE2(th1p1))';
M1p1_ip2 = inv(SE2(th1p1))*M{i+2}*inv(SE2(th1p1))';

M_mids1 = inv(SE2(th1s1))*M_mid_orgs1*inv(SE2(th1s1))';
M1s1_i = inv(SE2(th1s1))*M{i}*inv(SE2(th1s1))';
M1s1_is1 = inv(SE2(th1s1))*M{i-1}*inv(SE2(th1s1))';
M1s1_is2 = inv(SE2(th1s1))*M{i-2}*inv(SE2(th1s1))';

if isnan(M_mid)
    M_mid = eye(2,2);
end

y_org = sign(y1)*sqrt((U1(i,:) - Mid_pts)*((M1_i+M_mid)/2)*(U1(i,:) - Mid_pts)');
x_org1 = sign(x1)*sqrt((Mid_pts - U1(i-1,:))*((M1_is1+M_mid)/2)*(Mid_pts - U1(i-1,:))');
x_org2 = sign(x2)*sqrt((U1(i+1,:) - Mid_pts)*((M1_ip1+M_mid)/2)*(U1(i+1,:) - Mid_pts)');

y_orgp1 = sign(y1p1)*sqrt((U1p1(i+1,:) - Mid_ptsp1)*((M1p1_ip1+M_midp1)/2)*(U1p1(i+1,:) - Mid_ptsp1)');
x_org1p1 = sign(x1p1)*sqrt((Mid_ptsp1 - U1p1(i,:))*((M1p1_i+M_midp1)/2)*(Mid_ptsp1 - U1p1(i,:))');
x_org2p1 = sign(x2p1)*sqrt((U1p1(i+2,:) - Mid_ptsp1)*((M1p1_ip2+M_midp1)/2)*(U1p1(i+2,:) - Mid_ptsp1)');

y_orgs1 = sign(y1s1)*sqrt((U1s1(i-1,:) - Mid_ptss1)*((M1s1_is1+M_mids1)/2)*(U1s1(i-1,:) - Mid_ptss1)');
x_org1s1 = sign(x1s1)*sqrt((Mid_ptss1 - U1s1(i-2,:))*((M1s1_is2+M_mids1)/2)*(Mid_ptss1 - U1s1(i-2,:))');
x_org2s1 = sign(x2s1)*sqrt((U1s1(i,:) - Mid_ptss1)*((M1s1_i+M_mids1)/2)*(U1s1(i,:) - Mid_ptss1)');

delta_thtest = -(atan2(y_org,x_org1) + atan2(y_org,x_org2));
delta_thtestp1 = (atan2(y_orgp1,x_org1p1) + atan2(y_orgp1,x_org2p1));
delta_thtests1 = (atan2(y_orgs1,x_org1s1) + atan2(y_orgs1,x_org2s1));


% Now take derivative of this
if y_org ~= 0
    dy_org = (2*[0 1]*((M1_i+M_mid)/2)*(U1(i,:) - Mid_pts)'/(2*y_org));% + 0.5*(U1(i,:) - Mid_pts)*((dMy ))*(U1(i,:) - Mid_pts)'/(2*y_org);
    dx_org1 = 2*[1 0]*((M1_is1+M_mid)/2)*(Mid_pts - U1(i-1,:))'/(2*x_org1);
else
    dy_org = 0;
    dx_org1 = 0;
end

if y_orgp1 ~= 0
    dy_org_p1 = (2*[0 -0.5]*((M1p1_ip1+M_midp1)/2)*(U1p1(i+1,:) - Mid_ptsp1)'/(2*y_orgp1));% + 0.5*(U1p1(i+1,:) - Mid_ptsp1)*((dMyp1 ))*(U1p1(i+1,:) - Mid_ptsp1)'/(2*y_orgp1);
    dx_org1_p1 = 2*[-1 0]*((M1p1_i+M_midp1)/2)*(Mid_ptsp1 - U1p1(i,:))'/(2*x_org1p1);
else
    dy_org_p1 = 0;
    dx_org1_p1 = 0;
end

if y_orgs1 ~= 0
    dy_org_s1 = (2*[0 -0.5]*((M1s1_is1+M_mids1)/2)*(U1s1(i-1,:) - Mid_ptss1)'/(2*y_orgs1));% + 0.5*(U1s1(i-1,:) - Mid_ptss1)*((dMys1 ))*(U1s1(i-1,:) - Mid_ptss1)'/(2*y_orgs1);
    dx_org2_s1 = 2*[-1 0]*((M1s1_i+M_mids1)/2)*(U1s1(i,:) - Mid_ptss1)'/(2*x_org2s1);
else
    dy_org_s1 = 0;
    dx_org2_s1 = 0;
end


ddelta_thtest_y = ( dy_org*((1/x_org1 ) / (1 + (y_org/x_org1)^2)) + dy_org*(1/x_org2) / (1 + (y_org/x_org2)^2));
ddelta_thtest_yp1 = ( dy_org_p1*((1/x_org1p1 ) / (1 + (y_orgp1/x_org1p1)^2)) + dy_org_p1*(1/x_org2p1) / (1 + (y_orgp1/x_org2p1)^2));
ddelta_thtest_ys1 = ( dy_org_s1*((1/x_org1s1 ) / (1 + (y_orgs1/x_org1s1)^2)) + dy_org_s1*(1/x_org2s1) / (1 + (y_orgs1/x_org2s1)^2));

ddelta_thtest_x = ( -dx_org1*((y_org/(x_org1^2) ) / (1 + (y_org/x_org1)^2))) -dx_org1*((y_org/(x_org1^2) ) / (1 + (y_org/x_org1)^2));
ddelta_thtest_xp1 = ( -dx_org1_p1*((y_orgp1/(x_org1p1^2) ) / (1 + (y_orgp1/x_org1p1)^2)));
ddelta_thtest_xs1 = ( -dx_org2_s1*(y_orgs1/(x_org2s1^2)) / (1 + (y_orgs1/x_org2s1)^2)) ;
%%


% Parameterized Curvature (Note: the true curvature is after normalization)
kappa = delta_thtest/(1*(0.5*L_i+0.5*L_is1));

%-----Pretty sure Hossein said we'll only use Numerical method-----%
% Calculate the curvature from geodesic line (curvature calculated from the christoffel symbols)
% if strcmp(method, 'Analytical')
%     Gamma = Gamma1;
% else
    Gamma = Gamma2;
% end


k1_temp2 = ((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)*Gamma(:, :, 1)*((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)';
k2_temp2 = ((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)*Gamma(:, :, 2)*((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)';
% k3_temp2 = ((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)*Gamma(:, :, 3)*((r(i+1,:)-r(i,:))/2 + (r(i,:)-r(i-1,:))/2)';


k_temp = [k1_temp2 k2_temp2 ]*inv(T)*(Tr);

k1_temp = k_temp(1);
k2_temp = k_temp(2);

% Geodesic curvature which is how much a curve is away from being Geodesic
kappa_final2 =  (kappa - k2_temp/v_it);
kappa_final1 =  k1_temp/v_it;

% Calculate the gradient of length
shv = size(r,2);
for j = 1:shv
    
    dr = zeros(1,shv);
    dr(j) = 1;
    
    dL_i{j} = 1*(-dr*((M{i}+M{i+1} )/2)*(r(i+1,:)-r(i,:))' + (r(i+1,:)-r(i,:))*((M{i}+M{i+1} )/2)*(-dr') + 0.5*(r(i+1,:)-r(i,:))*((dM{i,j} ))*(r(i+1,:)-r(i,:))')/(2*L_i);
    dL_is1{j} = 1*(dr*((M{i-1}+M{i} )/2)*(r(i,:)-r(i-1,:))' + (r(i,:)-r(i-1,:))*((M{i-1}+M{i} )/2)*(dr') + 0.5*(r(i,:)-r(i-1,:))*(( dM{i,j}))*(r(i,:)-r(i-1,:))')/(2*L_is1);

end
    

dL_i = cell2mat(dL_i)*(Tr1);
dL_is1 = cell2mat(dL_is1)*(Tr1);

dL_ix = dL_i(1);
dL_iy = dL_i(2);

dL_is1x = dL_is1(1);
dL_is1y = dL_is1(2);



dkappa_y_temp = ((ddelta_thtest_yp1 + ddelta_thtest_ys1 + ddelta_thtest_y)*(1*(0.5*L_i+0.5*L_is1)) - delta_thtest*(dL_iy+dL_is1y)) /(1*(0.5*L_i+0.5*L_is1))^2 ;
dkappa_x_temp = (0*(ddelta_thtest_xp1 + ddelta_thtest_xs1 + ddelta_thtest_x)*(1*(0.5*L_i+0.5*L_is1)) - delta_thtest*(0.5*dL_ix+0.5*dL_is1x)) /(1*(0.5*L_i+0.5*L_is1))^2 ;


dkappa2 = [0 dkappa_y_temp]*inv(Tr1);
dkappa1 = [dkappa_x_temp 0]*inv(Tr1);


% Check the points on the rail in order to make the constraints

% Constraints check
% if ~isempty(cons)
%     % First need to check the node
%     for j = 1:length(cons.path_points)
%         
%         if ismember(i,cons.path_points{j})
% 
%             T_rail_normalized = cons.T_rail{j}/norm(cons.T_rail{j});
%             dkappa1 = dot(dkappa1,T_rail_normalized)*T_rail_normalized;
%             dkappa2 = dot(dkappa2,T_rail_normalized)*T_rail_normalized;
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

F = 0;



end
