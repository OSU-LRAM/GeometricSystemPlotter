function  [f,g1] = inertial_optimizer_calc(ytemp,s,P)
% [f,g1]

% Define the parameters
n = P.N;
% method = P.method; %-----Pretty sure Hossein said this isn't needed
y = [ytemp(1:n) ytemp(n+1:2*n)];
T_t = P.T;

%% Calculating cost and displacement per gait
p.phi = @(t) interp1( linspace(0,T_t,n+1), [y; y(1,:)], t);
for i=1:1:n-1

    velocityvalues(i,:)=n*(y(i+1,:)-y(i,:))/T_t;

end
velocityvalues(n,:)= n*(y(2,:)-y(1,:))/T_t;
p.dphi = @(t) interp1( linspace(0,T_t,n), [velocityvalues], t);

% Calculate the shape velocity
[net_disp_orig, net_disp_opt, PathLengthcost] = evaluate_displacement_and_cost1(s,p,[0, T_t],'interpolated','fixed_step',100);
lineint = net_disp_opt(1);


% Calculating the metric for each weaight point
P.M_t = s.metricfield.metric_eval.content.metric;
% P.M_t = s.metricfield.display.content.metric;  % Raw mwtric
alpha1 = s.grid.metric_eval{1};
alpha2 = s.grid.metric_eval{2};

P.alpha_field = {alpha1,alpha2};

% Find the metric
for i = 1:n
    M{i} = Granular_metric_calc2({y(i,1),y(i,2)},P.M_t,alpha1,alpha2);
end

% Define the metric
Metric = @(a,b) Granular_metric_calc2({a,b},P.M_t,alpha1,alpha2);

% Calculate the derivative of metric
shv = P.shv;
dM = {};
Gamma = {};
for i = 1:n
    [dM_temp,Gamma_temp] = metric_derivative_calc(y,i,P.M_t,P.alpha_field,shv,M{i});%-----Removed method from 3rd to last spot
    dM = [dM;dM_temp(i,:)];
    Gamma = [Gamma;Gamma_temp];
end

% s.DA_optimized{3} = s.height_optimized_corrected{3};
%% Height Function Calculation
for i=1:1:n
    
    height(i,3) = interpn(s.grid.eval{1,1}, s.grid.eval{2,1}, s.DA_optimized{1}, y(i,1), y(i,2));
    height(i,1) = 0;
    height(i,2) = 0;

end

% Calculate PathLength cost or Delta-V cos
%----- Pretty sure Hossein said that this switch isn't needed and to only
%use the Delta-v case-----%
% switch P.costType
%     case 'Delta-v'
        
        acc = zeros(n,1);
        for i = 3:n-2

            [grad_delta_v_temp,delta_v_temp,acc] = gradientstroke1(i,P,y,M,dM,[],shv,Gamma,acc,Metric);%-----removed method

            grad_delta_v(i,:) = grad_delta_v_temp';
            delta_v(i,:) = delta_v_temp;

        end

        % Point 1
        y1 = [[y(n-1,1) y(n-1,2)]; [y(n,1) y(n,2)];[y(1,1) y(1,2)];[y(2,1) y(2,2)];[y(3,1) y(3,2)]];
        M1 = [M(n-1); M(n);M(1);M(2);M(3)];
        dM1 = [dM(n-1,:);dM(n,:);dM(1,:);dM(2,:);dM(3,:)];
        Gamma1 = [Gamma(n-1); Gamma(n);Gamma(1);Gamma(2);Gamma(3)];
        [grad_delta_v_temp1,delta_v_temp1,acc1] = gradientstroke1(3,P,y1,M1,dM1,[],shv,Gamma1,acc,Metric); %-----Removed method
        grad_delta_v(1,:) = grad_delta_v_temp1';

        % Point 2
        y2 = [[y(n,1) y(n,2)]; [y(1,1) y(1,2)];[y(2,1) y(2,2)];[y(3,1) y(3,2)];[y(4,1) y(4,2)]];
        M2 = [M(n); M(1);M(2);M(3);M(4)];
        dM2 = [dM(n,:);dM(1,:);dM(2,:);dM(3,:);dM(4,:)];
        Gamma2 = [Gamma(n); Gamma(1);Gamma(2);Gamma(3);Gamma(4)];
        [grad_delta_v_temp2,delta_v_temp2,acc2] = gradientstroke1(3,P,y2,M2,dM2,[],shv,Gamma2,acc,Metric); %-----Removed method
        grad_delta_v(2,:) = grad_delta_v_temp2';

        % Point n
        yn = [[y(n-2,1) y(n-2,2)];[y(n-1,1) y(n-1,2)]; [y(n,1) y(n,2)];[y(1,1) y(1,2)];[y(2,1) y(2,2)]];
        Mn = [M(n-2);M(n-1); M(n);M(1);M(2)];
        dMn = [dM(n-2,:);dM(n-1,:);dM(n,:);dM(1,:);dM(2,:)];
        Gamman = [Gamma(n-2);Gamma(n-1); Gamma(n);Gamma(1);Gamma(2)];
        [grad_delta_v_tempn,delta_v_tempn,accn] = gradientstroke1(3,P,yn,Mn,dMn,[],shv,Gamman,acc,Metric); %-----Removed method
        grad_delta_v(n,:) = grad_delta_v_tempn;

        %Point n-1
        yns1 = [[y(n-3,1) y(n-3,2)];[y(n-2,1) y(n-2,2)]; [y(n-1,1) y(n-1,2)];[y(n,1) y(n,2)];[y(1,1) y(1,2)]];
        Mns1 = [M(n-3);M(n-2); M(n-1);M(n);M(1)];
        dMns1 = [dM(n-3,:);dM(n-2,:);dM(n-1,:);dM(n,:);dM(1,:)];
        Gammans1 = [Gamma(n-3);Gamma(n-2); Gamma(n-1);Gamma(n);Gamma(1)];
        [grad_delta_v_tempns1,delta_v_tempns1,accns1] = gradientstroke1(3,P,yns1,Mns1,dMns1,[],shv,Gammans1,acc,Metric); %-----Removed method
        grad_delta_v(n-1,:) = grad_delta_v_tempns1;

        % Total cost function (Delta-v)
        totaldelta_v = sum((delta_v)) + (delta_v_temp1) + (delta_v_temp2) + (delta_v_tempn) + (delta_v_tempns1);
        
        cost = totaldelta_v;
        grad = grad_delta_v;
        jacobianeqi = zeros(n,2);
%         [~,jacobianeqi] = pathlengthcost(y,M,dM,P);

%     case 'PathLength'
%         
%         [jacobianstroke,jacobianeqi] = pathlengthcost(y,M,dM,P);
%         
%         cost = PathLengthcost;
%         grad = jacobianstroke;
%         
% end
    
%% Jacobiandisp-jacobian for displacement produced by gait
for i=2:1:n-1
    jacobiandisp(i,:) = jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),height(i,:));
end
jacobiandisp(1,:) = jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),height(1,:));
jacobiandisp(n,:) = jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),height(n,:));

% Find the total Gradient
totaljacobian = jacobiandisp/cost-(((lineint)*grad)/(cost^2)) + 0.5*jacobianeqi;

dx = [totaljacobian(:)];

f = -lineint/cost;
g1 = -[totaljacobian(:)];

%-----Commenting out what looks like a bunch of plotting stuff-----%
% for i=1:n
%     G(i) = y(i,1);
%     H(i) = y(i,2);
%     W(i) = -totaljacobian(i,2);
%     V(i) = -totaljacobian(i,1);
% end
% 
% scale1=1;
% 
% clf(figure(2)) %%jacobianstroke
% figure(2)
% scale=0;
% quiver(G,H,V,W,'AutoScale','on')
% hold on
% plot(G,H,'-ro')
% hold on
% 
% axis equal
% hold off
% % pause(0.01)
% drawnow


end