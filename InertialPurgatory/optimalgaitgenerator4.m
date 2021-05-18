function y=optimalgaitgenerator4(s,dimension,npoints,a1,a2,a3,a4,lb,ub,direction,handles)
%%%%%%%%%%%%%%
% This function takes an input gait and runs fmincon to find the neareast locally 
% optimal gait

%Inputs:
%
%s: System file which contains the connection vector field, CCF's and
%   metric data
%dimension: Indicates the number of shape variables of the system
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
% a1: Values of the points forming the gait along the first shape dimension
% a2: Values of the points forming the gait along the second shape dimension
% lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% direction: Direction to optimize travel: 1=x,2=y,3=theta
% 
% 
% Outputs: 
% 
% y: Matrix whose values indicate coordinates of points which form the optimal gait
%%%%%%%%%%%%
n=npoints;

n=npoints;
P1(:,1)=a1(1,1:n)';
P1(:,2)=a2(1,1:n)';
P1(:,3)=a3(1,1:n)';
P1(:,4)=a4(1,1:n)';
%% Finding fourier coeffecients.
% The first step is to go from a direct transcription of the initial gait
% to a fourier based parametrization. 
% fa is a cell where the ith element contains the coefficients for the fourier parametrization along the ith direction 

t=1:1:npoints+1;
fa=cell(dimension);
% The bounds ensure the fourier series terms have the right period
options = fitoptions('fourier4');
options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi/n];
options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi/n];

for i=1:1:dimension
    fa{i}=fit(t',[P1(:,i);P1(1,i)],'fourier4',options);
end

%% The next step is to setup the fmincon call. 
% y0 is the marix of all fourier series coefficients that describe the
%   initial gait
% nonlcon is the function that imposes the constraints that all the points
%   stay inside the grid
% outfun is the function that updates the gait on the GUI after every iteration 

A=[];
b=[];
Aeq=[];
beq=[];

nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};%
 lb1=[];
 ub1=[];

for i=1:dimension
    for j=1:length(nu)
        y0(j,i)=fa{i}.(nu{j});
    end
end

 options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter','Algorithm','sqp','SpecifyObjectiveGradient',true,'CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,'MaxFunEvals',20000,'TolCon',10^-2,'OutputFcn', @outfun);
 [yf fval exitflag output]=fmincon(@(y) solvedifffmincon(y,s,n,dimension,direction,lb,ub),y0,A,b,Aeq,beq,lb1,ub1,@(y) nonlcon(y,s,n,dimension,lb,ub),options);


%% Getting point position values from the result of fmincon
% This section helps us go back to a direct transcription parametrization
% of the optimal gait from a fourier series parametrization. y is a column vector
% that contains coordinates of all points forming the optimized gait

for i=1:1:n
    for j=1:dimension
        y1(i,j)=yf(1,j)+yf(2,j)*cos(i*yf(end,j))+yf(3,j)*sin(i*yf(end,j))+yf(4,j)*cos(2*i*yf(end,j))+...
            +yf(5,j)*sin(2*i*yf(end,j))+yf(6,j)*cos(3*i*yf(end,j))+yf(7,j)*sin(3*i*yf(end,j))+...
            +yf(8,j)*cos(4*i*yf(end,j))+yf(9,j)*sin(4*i*yf(end,j));%+yf(10,j)*cos(5*i*yf(end,j))+yf(11,j)*sin(5*i*yf(end,j));%+yf(12,j)*cos(6*i*yf(end,j))+yf(13,j)*sin(6*i*yf(end,j));
    end    
end
y=y1(:);

%% Uncomment for plotting the optimized gait. Potentially useful while debugging.
% for i=1:n
%     xu(i)=y(i);
%     yu(i)=y(n+i);
%     zu(i)=y(2*n+i);
% end
% 
% for i=1:n`
%     xf(i)=P1(i,1);
%     yf(i)=P1(i,2);
% end
% 
% figure(11)
% hold on
% plot3(xu,yu,zu)

end

function [f,g]=solvedifffmincon(y,s,n,dimension,lb,ub)
%%%%%%%%%%%%%
% This function calculates efficiency (or displacement, if
% that is the objective function) and its gradient with respect to the coefficients obtained
% by the fourier series parametrization
% 
% Inputs: 
% 
% y: Matrix containing the Fourier series coefficients
% s: System file which contains the connection vector field, CCF's and
%   metric data
% dimension: Indicates the number of shape variables of the system
% n: The number of points desired in a direct transcription parametrization
%   of the gaits
% lb: Lower bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% 
% Outputs:
% 
% f: Objective function value (This is negative of efficiency by default, can be
%   changed to displacement)
% g: Gradient of the objective function
%%%%%%%%%%%%%

%% Obtaining points from fourier coefficients
% The first step is to obtain a direct transcription of the gait from the
% fourier series parametrization. y is the matrix of coordinate values of
% the points forming the gait.
afactor=0.001;
coeff=y;
for i=1:1:n
    for j=1:dimension
        y1(i,j)=y(1,j)+y(2,j)*cos(i*y(end,j))+y(3,j)*sin(i*y(end,j))+y(4,j)*cos(2*i*y(end,j))+...
            +y(5,j)*sin(2*i*y(end,j))+y(6,j)*cos(3*i*y(end,j))+y(7,j)*sin(3*i*y(end,j))+...
            +y(8,j)*cos(4*i*y(end,j))+y(9,j)*sin(4*i*y(end,j));%+y(10,j)*cos(5*i*y(end,j))+y(11,j)*sin(5*i*y(end,j));%+y(12,j)*cos(6*i*y(end,j))+y(13,j)*sin(6*i*y(end,j));%
    end    
end
clear y
y=y1;

%% Calculating cost and displacement per gait
velocityvalues=zeros(n,dimension);
g=1; %assign a time period for executing the gait
p.phi_def = @(t) interp1( linspace(0,g,n+1), [y; y(1,:)], t); % function parametrizing the gait as a function of time

for i=1:1:n-1
    velocityvalues(i,:)=n*(y(i+1,:)-y(i,:))/g; 
end
velocityvalues(n,:)=n*(y(1,:)-y(n,:))/g;
p.dphi_def = @(t) interp1( linspace(0,g,n), [velocityvalues], t); % Shape space velocity as a function of time


[net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, g],'interpolated','ODE'); % Call to the function that obtains displacement, cost and efficiency of a gait
lineint=net_disp_opt(direction); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 
totalstroke=cost; % Cost of executing the gait ones

% If efficiency is negative, reversing the order of points so that
% efficiency is positive
if lineint<0
    lineint=-lineint;
    ytemp=y;
    for i=1:n
        y(i,:)=ytemp(n+1-i,:);
    end
    invert=1;
else
    lineint=lineint;
    invert=0;
end

%% Preliminaries for gradient calculation
% Preallocating memory for variables which we will need in further
% calculation 
yvalues=cell(n,dimension); % Cell representation of the coordinates of all points forming the gait
interpstateccf=cell(1,dimension); % Variable which will store the ccf function grid used for interpolation
interpmetricgrid=cell(1,dimension);  % Variable which will store the metric grid used for interpolation
ccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store ccf function at each point
metric1=zeros(n,dimension,dimension);% Variable which will store metric at each point in the form of a matrix
metric = repmat({zeros(dimension)},[n 1]); % Variable which stores the metric at each point in the form of a 2x2 matrix
metricgrad1=zeros(n,dimension,dimension,dimension);% Variable which will store gradient of metric at each point in the form of a matrix
metricgrad = repmat({zeros(dimension)},[n,dimension]);% Variable which will store gradient of metric at each point in the form of a matrix

% Interpolation to calculate all the variables needed for gradient
% calculation
for i=1:1:n
    for j=1:1:dimension
        yvalues{i,j}=y(i,j);
    end

    for j=1:1:dimension
        interpstateccf{j}=s.grid.eval{j,1};
        interpmetricgrid{j}=s.grid.metric_eval{j,1};
    end
end


for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y(:,1),y(:,2),y(:,3),y(:,4),'spline');
end


for j=1:1:dimension
    for k=1:1:dimension
        metric1(:,j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y(:,1),y(:,2),y(:,3),y(:,4),'spline');
    end
end
% metric{i}=eye(2);
for i=1:n
    for j=1:1:dimension
       for k=1:1:dimension
           metric{i}(j,k)=metric1(i,j,k);
       end
    end
end

for l=1:1:dimension
    for m=1:1:dimension
        if m==l
           y2(:,m)=y(:,m)+afactor*ones(length(y),1);
           y1(:,m)=y(:,m)-afactor*ones(length(y),1);
        else
           y2(:,m)=y(:,m);
           y1(:,m)=y(:,m);
        end
    end
    for j=1:1:dimension
        for k=1:1:dimension
            metricgrad1(:,l,j,k)=(interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y2(:,1),y2(:,2),y(:,3),y(:,4),'spline')-interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y1(:,1),y1(:,2),y(:,3),y(:,4),'spline'))/(2*afactor);
        end
    end
    for i=1:n
        for j=1:1:dimension
            for k=1:1:dimension
                metricgrad{i,l}(j,k)=metricgrad1(i,l,j,k);
            end
        end
    end
end

%l is the vector containing metric weighted distances between neighbouring
%points

for i=1:1:n-1
    l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
end
l(n)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');

%% Jacobianstroke is the gradient of cost. 
%Contrigrad is the contribution to the gradient due to the metric changing

jacobianstroke = zeros(n,dimension);
contrigrad=zeros(n,dimension);

for i=1:n-1
    delp{i}=y(i+1,:)-y(i,:); % delp{i} is the vector joining the (i+1)th point to the ith point 
end
delp{n}=y(1,:)-y(n,:);

for i=2:n-1
    for j=1:dimension
        %Contrigrad is the contribution to the gradient due to the metric changing
        contrigrad(i,j)=0.5*delp{i}*metricgrad{i,j}*delp{i}'/(2*l(i))+0.5*delp{i-1}*metricgrad{i,j}*delp{i-1}'/(2*l(i-1)); 
    end
    % Total gradient is the result of distance changing due to movement of point and the metric changing due to movement of the point
    jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+...
        +((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+contrigrad(i,:); 
end

% Calculation for the 1st point and last point have to be done outside the
% loop as the (i+1)th point for the last point is the first point and
% (i-1)th point for the first point is the last point
for j=1:dimension
    contrigrad(1,j)=0.5*delp{1}*metricgrad{1,j}*delp{1}'/(2*l(1))+0.5*delp{n}*metricgrad{1,j}*delp{n}'/(2*l(n));
end
jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+...
    +((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+contrigrad(1,:);

for j=1:dimension
    contrigrad(n,j)=0.5*delp{n}*metricgrad{n,j}*delp{n}'/(2*l(n))+0.5*delp{n-1}*metricgrad{n,j}*delp{n-1}'/(2*l(n-1));
end
jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+...
    +((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+contrigrad(n,:);   


%% Jacobiandisp is the gradient of displacement.
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. It's input arguments are the coordinates of 
% the (i-1)th, ith and (i+1)th point, CCF value at point i and the dimension of     
% the shape space (dimension)

jacobiandisp = zeros(n,dimension);
for i=2:1:n-1
    jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
end
jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);

%% Jacobianeqi is the concentration gradient. 
%It is the term that keeps points eqi distant from each other and prevents crossover of gait.

jacobianeqi = zeros(n,dimension);
for i=2:n-1;
    len=sqrt((y(i+1,:)-y(i-1,:))*((metric{i-1}+metric{i+1})/2)*(y(i+1,:)-y(i-1,:))'); % metric weighted length between point (i-1) and (i+1)
    midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2))/2; % location of midpoint of the line segment joining point (i-1) and (i+1)
    betacos=(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i})/2))'/(l(i-1)*len);
    xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*l(i-1)*betacos/len; %projection of ith point onto the line joining the (i-1)th and (i+1)th points
    jacobianeqi(i,:)=midpoint-xhat; % gradient of the ith point is equal to the difference between the midpoint and the projection of ith point
end

    len=sqrt((y(2,:)-y(n,:))*((metric{2}+metric{n})/2)*(y(2,:)-y(n,:))');
    midpoint=y(n,:)+((y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2))/2;
    betacos=(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*((y(1,:)-y(n,:))*sqrtm((metric{n}+metric{1})/2))'/(l(n)*len);
    xhat=y(n,:)+(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*l(n)*betacos/len;
    jacobianeqi(1,:)=midpoint-xhat;

    len=sqrt((y(1,:)-y(n-1,:))*((metric{1}+metric{n-1})/2)*(y(1,:)-y(n-1,:))');
    midpoint=y(n-1,:)+((y(1,:)-y(n-1,:))*sqrtm((metric{1}+metric{n-1})/2))/2;
    betacos=(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*((y(n,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{n})/2))'/(l(n-1)*len);
    xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*l(n-1)*betacos/len;
    jacobianeqi(n,:)=midpoint-xhat;

%% changey/dcoeff tells us how much each point moves when a fourier series variable is changed
% chy is a cell with as many entries as the dimension of the shape space
% ith element of chy is a matrix where the (j,k)th entry tells us the change in the ith coordinate
% of the kth point of the gait resulting from a unit change in the jth
% fourier coefficient corresponding to the ith dimension of the shape space

chy=cell(dimension,1);
for i=1:1:dimension
    for j=1:1:n
        chy{i}(:,j)=[1;cos(j*coeff(end,i));sin(j*coeff(end,i));cos(2*j*coeff(end,i));sin(2*j*coeff(end,i));cos(3*j*coeff(end,i));sin(3*j*coeff(end,i));cos(4*j*coeff(end,i));sin(4*j*coeff(end,i))];%cos(5*j*coeff(end,i));sin(5*j*coeff(end,i))];%;cos(6*j*coeff(end,i));sin(6*j*coeff(end,i))];%
    end
end



%% properly ordering gradients depending on wether lineint was negative or positive
if invert==0
    jacobiandisptemp=jacobiandisp;
    jacobianstroketemp=jacobianstroke;
    jacobianeqitemp=jacobianeqi;
    jacobiandisp=jacobiandisp;
    jacobianstroke=jacobianstroke;
    jacobianeqi=jacobianeqi;
else
        jacobiandisptemp=jacobiandisp;
        jacobianstroketemp=jacobianstroke;
        jacobianeqitemp=jacobianeqi;
    for i=1:1:n
        jacobiandisp(i,:)=jacobiandisptemp(n+1-i,:);
        jacobianstroke(i,:)=jacobianstroketemp(n+1-i,:);
        jacobianeqi(i,:)=jacobianeqitemp(n+1-i,:);
    
    end
end

%% fourier series version of all gradients

% The line below calculates the total gradient in a direct transcription
% parametrization
totaljacobian=jacobiandisp/totalstroke-lineint*jacobianstroke/totalstroke^2+jacobianeqi;

% We then obtain gradients in a fourier series parametrization by
% projecting the gradients from the direct transcription space onto the
% fourier coefficient space
for i=1:1:dimension
    for j=1:1:9 
        jacobiandispfourier(j,i)=chy{i}(j,:)*jacobiandisp(:,i);
        jacobianstrokefourier(j,i)=chy{i}(j,:)*jacobianstroke(:,i);
        jacobianeqifourier(j,i)=chy{i}(j,:)*jacobianeqi(:,i);
        totaljacobianfourier(j,i)=chy{i}(j,:)*totaljacobian(:,i);
    end
end
    





%% Variables useful while debugging for flaws in gradient calculations

% for i=1:n
%     for j=1:1:dimension
%         totaljacobianc(i,j)=chy{j}(:,i)'*totaljacobianfourier(:,j);
%         jacobiandispc(i,j)=chy{j}(:,i)'*jacobiandispfourier(:,j);
%         jacobianstrokec(i,j)=chy{j}(:,i)'*jacobianstrokefourier(:,j);
%     end
% end
%         totaljacobianctemp=totaljacobianc;
%         jacobiandispctemp=jacobiandispc;
%         jacobianstrokectemp=jacobianstrokec;
% for i=1:1:n
%     jacobiandispc(i,:)=jacobiandispctemp(n+1-i,:);
%     jacobianstrokec(i,:)=jacobianstrokectemp(n+1-i,:);
%     totaljacobianc(i,:)=totaljacobianctemp(n+1-i,:);
% end

%% minimizing negative of efficiency(or displacement)
 f=-lineint/(totalstroke);
% f=-lineint;
if nargout>1
%     g=-totaljacobian(:);
    g=[-totaljacobianfourier;zeros(1,dimension)];
end

%% Debugging and plotting
% This section was written up for the sole purpose of helping with
% debugging flaws in the optimizer. Appropriate sections can be uncommented
% to plot how the different gradients look during the optimization process.

% for i=1:n
%     G(i)=y(i,1);
%     H(i)=y(i,2);
%     P(i)=y(i,3);
%     I(i)=1*jacobianstroketemp(i,1);
%     J(i)=1*jacobianstroketemp(i,2);
%     I(i)=1*jacobianstroke(i,1);
%     J(i)=1*jacobianstroke(i,2);
%     N(i)=1*jacobianstroke(i,3);
%     K(i)=jacobiandisptemp(i,1);
%     K(i)=jacobiandisp(i,1);
%     K1(i) = jacobiandispc(i,1);
%     K1(i)=jacobianforward(i,1);
%     L(i)=jacobiandisptemp(i,2);
%     L(i)=jacobiandisp(i,2);
%     L1(i)=jacobiandispc(i,2);
%     L1(i)=jacobianforward(i,2);
%     Q(i)=jacobiandisp(i,3);
%     B(i)=totaljacobianc(i,1);
%     B1(i)=totaljacobian(i,1);
%     C(i)=totaljacobianc(i,2);
%     C1(i)=totaljacobian(i,2);
%     D(i)=totaljacobian(i,3);
%     O(i)=jacobianeqi(i,1);
%     S(i)=jacobianeqi(i,2);
%     R(i)=jacobianeqi(i,3);
%     ccfx(i)=ccf(i,1);
%     ccfy(i)=ccf(i,2); 
%     ccfz(i)=ccf(i,3);
%     ccfcrosscheck(i)=ccf(i,:)*jacobiandisp(i,:)';
% end
% % % % 
% % % % 
% clf(figure(6)) %%jacobiandisp
% figure(6)
% scale=0;
% scale1=1;
% quiver(G,H,K,L,scale1)
% hold on
% quiver(G,H,K1,L1,scale1)
% % quiver(G,H,ccfx,ccfy,scale1)
% plot(G,H)
% axis equal
% hold off
% 
% % clf(figure(9)) %%jacobianforward
% % figure(9)
% % scale=0;
% % scale1=1;
% % quiver(G,H,K1,L1,scale1)
% % hold on
% % % quiver(G,H,P,ccfx,ccfy,ccfz,scale1)
% % plot(G,H)
% % axis equal
% % hold off
% % 
% 
% 
% % % 
% clf(figure(7)) %%jacobianstroke
% figure(7)
% scale=0;
% quiver(G,H,K,L,scale)
% quiver(G,H,I,J,scale1)
% hold on
% plot(G,H)
% axis equal
% hold off
% % % % 
% clf(figure(8)) %%% totaljacobian
% figure(8)
% scale=0;
%  quiver(G,H,B,C,scale1)
% hold on
% quiver(G,H,B1,C1,scale1)
% 
% plot(G,H)
% axis equal
% hold off
% % 
% % clf(figure(3)) %%%jacobianeqi
% % figure(3)
% % scale1=1;
% % quiver3(G,H,P,O,S,R,scale1)
% % hold on
% % plot3(G,H,P)
% % axis equal
% % hold off
% % 
% % 
% % % figure(5)
% % % scale=0;
% % % quiver(G,H,B1,C1,scale1)
% % % hold on
% % % plot(G,H)
% % % axis equal
% % % hold off
% % min(abs(totaljacobian))
% % max(abs(totaljacobian))
% % pause(0.1)
% 
% % 
% % for i=1:1:length(ccf)
% %     ccfcheckper(i)=(ccfcheck(i)-ccf(i))/ccf(i);
% % end
% 
% figure(5)
% plot(y(:,1),y(:,2))
% axis equal
% 
%  pause(0.1)



end

function a=jacobiandispcalculator3(p1,p2,p3,ccf,dimension)
%%%%%%%%%
%
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. 
% It's input arguments are the coordinates of the (i-1)th, ith and (i+1)th point,     
% CCF value at point i(ccf) and the dimension of the shape space (dimension)
%
%%%%%%%%%

l1=0; % variable for calculating length of the line segment joining the (i-1)th point with the (i+1)th point
for i=1:1:dimension
    l1=l1+(p1(i)-p3(i))^2;
    base(1,i)=p3(i)-p1(i); % vector connecting the (i-1)th point and (i+1)th point  
end
l=sqrt(l1); % length of the line segment joining the (i-1)th point with the (i+1)th point

for i=1:1:dimension
    jacobian(1,i)=0;
    perp1=zeros(1,dimension);
    perp1(i)=1;
    %parcomp=base*perp1'/norm(base);
    %perp1-parcomp*base/norm(base);  %%recheck again
    perp=perp1;% Unit vector along the ith direction
    % The for loop below calculates the gradient along the ith direction by
    % treating the CCF as 2 forms. A specific (j,k) represents a component of the 2 form 
    for j=1:dimension-1 
        for k=1:dimension-j
            veca=zeros(1,dimension);
            vecb=zeros(1,dimension);
            veca(j)=1;
            vecb(j+k)=1;
            f=(j-1)*dimension-(j*(j-1))/2+k;
            jacobian(1,i)=jacobian(1,i)+0.5*ccf(f)*((veca*perp')*(vecb*base')-(vecb*perp')*(veca*base'));
        end
    end
end

a=jacobian;

end

function [A,Aeq]=nonlcon(y,s,n,dimension,lb,ub)
%%%%%%%%% 
%
%This function imposes the nonlinear constraint that all the points forming the gait stay within bounds
%
%Inputs:
%
%y: Fourier series coefficients that describe the gait
%s: System file which contains the connection vector field, CCF's and
%   metric data
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
%dimension: Indicates the number of shape variables of the system
%lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
%ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% 
%%%%%%%%%

% The first step is to obtain a direct transciption parametrization of the gait from the 
% fourier series parametrization
for i=1:1:n+1
    for j=1:dimension
        y1(i,j)=y(1,j)+y(2,j)*cos(i*y(end,j))+y(3,j)*sin(i*y(end,j))+y(4,j)*cos(2*i*y(end,j))+...
            +y(5,j)*sin(2*i*y(end,j))+y(6,j)*cos(3*i*y(end,j))+y(7,j)*sin(3*i*y(end,j))+...
            +y(8,j)*cos(4*i*y(end,j))+y(9,j)*sin(4*i*y(end,j));%+y(10,j)*cos(5*i*y(end,j))+y(11,j)*sin(5*i*y(end,j));%+y(12,j)*cos(6*i*y(end,j))+y(13,j)*sin(6*i*y(end,j));
    end    
end

y2=y1(:);

b=length(y2);


% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A=[A1;A2];

% Aeq ensures that fmincon doesn't alter the period of fourier series
% components
Aeq=y(end,:)-2*pi/n*ones(size(y(end,:)));

end

function stop=outfun(y,optimValues,state,handles)
%%%%%%%%% 
%
%This function plots the current state of the gait on the sysplotter GUI
%after every iteration of the optimizer
%
%%%%%%%%% 
n=100;
dimension=length(y(1,:));

for thisAxes = [1:numel(handles.plot_thumbnails.Children)]
    
    axes(handles.plot_thumbnails.Children(thisAxes));
    % The if else statement below deletes gaits 2 iterations after they have been plotted
    if optimValues.iteration>2
        children=get(gca,'children');
        delete(children(2));
    else
    end

    % The if else statement below fades the gait plotted during the previous iteration
    if optimValues.iteration>1
        children=get(gca,'children');
        children(1).Color=[0.5 0.5 0.5];
        children(1).LineWidth=4;
    else
    end
    % 
    % % The if else statement below plots the gait after every iteration
     if optimValues.iteration>0
        for i=1:1:n+1
            for j=1:dimension
                y1(i,j)=y(1,j)+y(2,j)*cos(i*y(end,j))+y(3,j)*sin(i*y(end,j))+y(4,j)*cos(2*i*y(end,j))+...
                    +y(5,j)*sin(2*i*y(end,j))+y(6,j)*cos(3*i*y(end,j))+y(7,j)*sin(3*i*y(end,j))+...
                    +y(8,j)*cos(4*i*y(end,j))+y(9,j)*sin(4*i*y(end,j));%+y(10,j)*cos(5*i*y(end,j))+y(11,j)*sin(5*i*y(end,j));%+y(12,j)*cos(6*i*y(end,j))+y(13,j)*sin(6*i*y(end,j));
            end    
        end
        hold on
        handle1=plot3(y1(:,1),y1(:,2),y1(:,3),'k','linewidth',3);
    else
    end
end


pause(0.1)
stop=false;
end

function [net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,tspan,ConnectionEval,IntegrationMethod,resolution)
% Evaluate the displacement and cost for the gait specified in the
% structure GAIT when executed by the system described in the structure
% S.
%
% S should be a sysplotter 's' structure loaded from a file
% sysf_FILENAME_calc.mat (note the _calc suffix)
%
% P should be a structure with fields "phi_def" and "dphi_def", returning a
% vector of shapes and shape velocities respectively. If it is not
% convenient to analytically describe the shape velocity function,
% gait.dphi should be defined as 
%
% p.dphi =  @(t) jacobianest(@(T) p.phi (T),t)
%
% as is done automatically by sysplotter, but note that this will be slower
% than specifying a function directly
%
% ConnectionEval can specify whether the local connection should be generated from
% its original function definiton, or by interpolation into the evaluated
% matrix, but is optional. Valid options are 'functional' or
% 'interpolated'. Defaults to "interpolated", which significantly faster
% when calculating the local connection or metric from scratch takes
% appreciable computational time
%
% IntegrationMethod can specify whether ODE45 or a fixed point
% (euler-exponential) integration method should be employed. Defaults to
% ODE, fixed point code is experimental.
%
% RESOLUTION specifies the number of points for fixed-point resolution
% evaluation. A future option may support autoconvergence, but ODE
% performance with interpolated evaluation appears to be fast enough that
% refinement of fixed-point performance is on hold.
	

	% if no ConnectionEval method is specified, default to interpolated
	if ~exist('ConnectionEval','var')
		ConnectionEval = 'interpolated';
	end
    
    % if no IntegrationMethod is specified, default to ODE
	if ~exist('IntegrationMethod','var')
		IntegrationMethod = 'ODE';
	end

    % if no resolution is specified, default to 100 (this only affects
    % fixed_step integration)
	if ~exist('resolution','var')
		resolution = 100;
	end

    
    
	switch IntegrationMethod
		
		case 'fixed_step'
			
			[net_disp_orig, cost] = fixed_step_integrator(s,p,tspan,ConnectionEval,resolution);
        
        case 'ODE'

            % Calculate the system motion over the gait
            sol = ode45(@(t,y) helper_function(t,y,s,p,ConnectionEval),tspan,[0 0 0 0]');

            % Extract the final motion
            disp_and_cost = deval(sol,tspan(end));

            net_disp_orig = disp_and_cost(1:3);
            cost = disp_and_cost(end);

            % Convert the final motion into its representation in optimal
            % coordinates
            startshape = p.phi_def(0);
            startshapelist = num2cell(startshape);
            beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'spline');
            net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
                -sin(beta_theta) cos(beta_theta) 0;...
                0 0 1]*net_disp_orig;
            
        otherwise
			error('Unknown method for integrating motion');
	end

	
	% Convert the final motion into its representation in optimal
	% coordinates
	startshape = p.phi_def(0);
	startshapelist = num2cell(startshape);
	beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'spline');
	net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
		-sin(beta_theta) cos(beta_theta) 0;...
		0 0 1]*net_disp_orig;

	
end

% Evaluate the body velocity and cost velocity (according to system metric)
% at a given time
function [xi, dcost] = get_velocities(t,s,gait,ConnectionEval)

	% Get the shape and shape derivative at the current time
	shape = gait.phi_def(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi_def(t);
	
	
	% Get the local connection and metric at the current time, in the new coordinates	
	switch ConnectionEval
		case 'functional'
			
			A = s.A_num(shapelist{:})./s.A_den(shapelist{:});
			
			M = s.metric(shapelist{:});

		case 'interpolated'
			
			A = -cellfun(@(C) interpn(s.grid.eval{:},C,shapelist{:},'spline'),s.vecfield.eval.content.Avec);
			
			M =  cellfun(@(C) interpn(s.grid.metric_eval{:},C,shapelist{:},'spline'),s.metricfield.metric_eval.content.metric);
			
		otherwise
			error('Unknown method for evaluating local connection');
	end
	
	% Get the body velocity at the current time
	t;
    xi = - A * dshape(:);

	% get the cost velocity
	dcost = sqrt(dshape(:)' * M * dshape(:));
	
end


% Function to integrate up system velocities using a fixed-step method
function [net_disp_orig, cost] = fixed_step_integrator(s,gait,tspan,ConnectionEval,resolution)

	% Duplicate 'resolution' to 'res' if it is a number, or place res at a
	% starting resolution if an automatic convergence method is selected
	% (automatic convergence not yet enabled)
	default_res = 100;
	if isnumeric(resolution)
		res = resolution;
	elseif isstr(resolution) && strcmp(resolution,'autoconverge')
		res = default_res;
	else
		error('Unexpected value for resolution');
	end
	
	% Generate the fixed points from the time span and resolution
	tpoints = linspace(tspan(1),tspan(2),res);
	tsteps = gradient(tpoints);

	% Evaluate the velocity function at each time
	[xi, dcost] = arrayfun(@(t) get_velocities(t,s,gait,ConnectionEval),tpoints,'UniformOutput',false);
	
	
	%%%%%%%
	% Integrate cost and displacement into final values
	
	%%
	% Exponential integration for body velocity
	
	% Exponentiate each velocity over the corresponding time step
	expXi = cellfun(@(xi,timestep) se2exp(xi*timestep),xi,num2cell(tsteps),'UniformOutput',false);
	
	% Start off with zero position and displacement
	net_disp_matrix = eye(size(expXi{1}));
	
	% Loop over all the time steps from 1 to n-1, multiplying the
	% transformation into the current displacement
	for i = 1:(length(expXi)-1)
		
		net_disp_matrix = net_disp_matrix * expXi{i};
		
	end
	
	% De-matrixafy the result
	g_theta = atan2(net_disp_matrix(2,1),net_disp_matrix(1,1));
	g_xy = net_disp_matrix(1:2,3);
	
	net_disp_orig = [g_xy;g_theta];
	
	%%
	% Trapezoidal integration for cost
	dcost = [dcost{:}];
	cost = trapz(tpoints,dcost);

end


% Function to evaluate velocity and differential cost at each time for ODE
% solver
function dX = helper_function(t,X,s,gait,ConnectionEval)

	% X is the accrued displacement and cost

	[xi, dcost] = get_velocities(t,s,gait,ConnectionEval);
		
	% Rotate body velocity into world frame
	theta = X(3);
	v = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
		
	% Combine the output
	dX = [v;dcost];
	

end

function expXi = se2exp(xi)

	% Make sure xi is a column
	xi = xi(:);

	% Special case non-rotating motion
	if xi(3) == 0
		
		expXi = [eye(2) xi(1:2); 0 0 1];
		
	else
		
		z_theta = xi(3);
		
		z_xy = 1/z_theta * [sin(z_theta), 1-cos(z_theta); cos(z_theta)-1, sin(z_theta)] * xi(1:2);
		
		expXi = [ [cos(z_theta), -sin(z_theta); sin(z_theta), cos(z_theta)], z_xy;
			0 0 1];
		
	end


end

function [g_end_orig,g_end_opt, cost_end] = extract_displacement_and_cost(datafile)
% Extract the displacement and cost data from a sysf_...shchf_....mat file

% Load the target file
load(datafile,'p')

% Prime arrays to hold the net displacement (in original and optimal
% coordinates) and cost from each shape change in the file. p.G_locus_full is
% single-level cell array of structures, each of which holds the
% information for one gait (with all segments concatenated)
g_end_orig = zeros(numel(p.G_locus_full),3);
g_end_opt = g_end_orig;
cost_end = zeros(numel(p.G_locus_full,1)); % If distance metric was not specified, euclidean metric in the parameters was assumed

% Loop over each shape change
for i = 1:numel(p.G_locus_full)
	
	% Extract the final values for the relevant parameters
	g_end_orig(i,:) = p.G_locus_full{i}.G(end,:); 
	g_end_opt(i,:) = p.G_locus_full{i}.G_opt(end,:); 
	cost_end(i) = p.G_locus_full{i}.S(end);
end
end