function y=optimalgaitgenerator(s,dimension,npoints,a1,a2,lb,ub)

% load('sysf_honey_swimmer_calc.mat');
% npoints=100;
% dimension=2;
% y=optimalgaitgenerator(s,2,100)


n=npoints;

% for i=1:1:n
%     P1(i,1)=1*cos((i-1)*2*pi/n);
%     P1(i,2)=1*sin((i-1)*2*pi/n);
%     for j=3:1:dimension
%         P1(i,j)=0;
%     end
% end
P1(:,1)=a1(1,1:n)';
P1(:,2)=a2(1,1:n)';

% for i=1:1:n/2
%     P1(i,1)=-0.5+0.5*cos((i-1)*4*pi/n);
%     P1(i,2)=+0.5*sin((i-1)*4*pi/n);
%     for j=3:1:dimension
%         P1(i,j)=0.0;
%     end
% end
% 
% for i=n/2+1:1:n
%     P1(i,1)=0.5-0.5*cos(-(i-1)*4*pi/n);
%     P1(i,2)=-0.5*sin(-(i-1)*4*pi/n);
%     for j=3:1:dimension
%         P1(i,j)=0.0;
%     end
%end

P0=P1(:);
A=[];
b=[];
Aeq=[];
beq=[];
% lb=-2.2*ones(200,1);
% ub=2.2*ones(200,1);
nonlcon=[];
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter','Algorithm','active-set','GradObj','on','TolX',10^-4,'TolFun',10^-6,'MaxIter',4000,'MaxFunEvals',20000);
 [y fval exitflag output]=fmincon(@(y) solvedifffmincon(y,s,n,dimension),P0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%  [y fval exitflag output]=fmincon(@(y) solvedifffmincon1(y,s,n),P0,A,b,Aeq,beq,lb,ub,nonlcon,options);

for i=1:n
    xf(i)=y(i);
    yf(i)=y(n+i);
end

% for i=1:n
%     xf(i)=P1(i,1);
%     yf(i)=P1(i,2);
% end

% figure(10)
% hold on
% plot(xf,yf,'r')
end

function [f,g]=solvedifffmincon(y,s,n,dimension)
afactor=0.001;
y1=y;
clear y
y=zeros(n,dimension);
for j=1:1:dimension
    y(:,j)=y1((j-1)*n+1:n*j);
end

% pointvalues=[y1(1:n)*abar,y1(n+1:2*n)*bbar,y1(2*n+1:3*n)*cbar];
pointvalues=y;

%% Calculating cost and displacement per gait
g=10;
p.phi = @(t) interp1( linspace(0,g,n+1), [pointvalues; pointvalues(1,:)], t);
velocityvalues = zeros(n-1,size(pointvalues,2));
for i=1:1:n-1
    velocityvalues(i,:)=n*(pointvalues(i+1,:)-pointvalues(i,:))/g;
end
velocityvalues(n,:)=n*(pointvalues(1,:)-pointvalues(n,:))/g;
p.dphi = @(t) interp1( linspace(0,g,n), [velocityvalues], t);


[net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, g],'interpolated','fixed_step',100);
lineint=net_disp_opt(1);
totalstroke=cost;

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


%% Preliminaries for calculation

height = zeros(n,dimension*(dimension-1)/2);
metric = repmat({zeros(dimension)},[n 1]);
yvalues1 = cell(dimension,1);
yvalues2 = cell(dimension,1);
metricgrad = repmat({zeros(dimension)},[n,dimension]);

for i=1:1:n
%     for j=1:1:dimension
%         yvalues{j}=y(i,j);
%     end
    
    yvalues = num2cell(y(i,:));

%     for j=1:1:dimension
%         interpstateheight{j}=s.grid.eval{j,1};
%         interpmetricgrid{j}=s.grid.metric_eval{j,1};
%     end

    interpstateheight = s.grid.eval(:,1);
    interpmetricgrid = s.grid.metric_eval(:,1);
    
    for j=1:size(height,2)
        height(i,j)=interpn(interpstateheight{:},s.DA_optimized{1,j},yvalues{:},'cubic');
    end
    
%    metricsize=dimension*dimension;
    
    for j=1:1:dimension
        for k=1:1:dimension
            metric{i}(j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},yvalues{:},'cubic');
        end
    end
% metric{i}=eye(2);
    
    for l=1:1:dimension
        for m=1:1:dimension
            if m==l
                yvalues2{m}=y(i,m)+afactor;
                yvalues1{m}=y(i,m)-afactor;
            else
                yvalues2{m}=y(i,m);
                yvalues1{m}=y(i,m);
            end
        end
            for j=1:1:dimension
                for k=1:1:dimension
                    metricgrad{i,l}(j,k)=(interpn(interpmetricgrid{:}, s.metricfield.metric_eval.content.metric{j,k}, yvalues2{:},'cubic')-interpn(interpmetricgrid{:}, s.metricfield.metric_eval.content.metric{j,k}, yvalues1{:},'cubic'))/(2*afactor);
                end
            end
        
    end
end

% for i=2:1:n-1
%     lengthjacdisp(i)= sqrt((y(i-1,1)-y(i+1,1))^2+(y(i-1,2)-y(i+1,2))^2+(y(i-1,3)-y(i+1,3))^2);
% end
% lengthjacdisp(1)=sqrt((y(n,1)-y(2,1))^2+(y(n,2)-y(2,2))^2+(y(n,3)-y(2,3))^2);
% lengthjacdisp(n)=sqrt((y(n-1,1)-y(1,1))^2+(y(n-1,2)-y(1,2))^2+(y(n-1,3)-y(1,3))^2);

pointvalues=y;

for i=1:1:n-1
    l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
end
l(n)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');

%% Jacobianstroke

jacobianstroke = zeros(n,dimension);
contrigrad=zeros(n,dimension);

for i=1:n-1
    delp{i}=y(i+1,:)-y(i,:);
end
delp{n}=y(1,:)-y(n,:);

for i=2:n-1
    for j=1:dimension
        contrigrad(i,j)=0.5*delp{i}*metricgrad{i,j}*delp{i}'/(2*l(i))+0.5*delp{i-1}*metricgrad{i,j}*delp{i-1}'/(2*l(i-1));
    end
    jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+contrigrad(i,:);%+0.5*[delp{i}*metricgradx{i}*delp{i}',delp{i}*metricgrady{i}*delp{i}',delp{i}*metricgradz{i}*delp{i}']/(2*l(i))+0.5*[delp{i-1}*metricgradx{i}*delp{i-1}',delp{i-1}*metricgrady{i}*delp{i-1}',delp{i-1}*metricgradz{i}*delp{i-1}']/(2*l(i-1));
end

for j=1:dimension
    contrigrad(1,j)=0.5*delp{1}*metricgrad{1,j}*delp{1}'/(2*l(1))+0.5*delp{n}*metricgrad{1,j}*delp{n}'/(2*l(n));
end
jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+contrigrad(1,:);%+0.5*[delp{1}*metricgradx{1}*delp{1}',delp{1}*metricgrady{1}*delp{1}',delp{1}*metricgradz{1}*delp{1}']/(2*l(1))+0.5*[delp{n}*metricgradx{1}*delp{n}',delp{n}*metricgrady{1}*delp{n}',delp{n}*metricgradz{1}*delp{n}']/(2*l(n));   

for j=1:dimension
    contrigrad(n,j)=0.5*delp{n}*metricgrad{n,j}*delp{n}'/(2*l(n))+0.5*delp{n-1}*metricgrad{n,j}*delp{n-1}'/(2*l(n-1));
end
jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+contrigrad(n,:);%+0.5*[delp{n}*metricgradx{n}*delp{n}',delp{n}*metricgrady{n}*delp{n}',delp{n}*metricgradz{n}*delp{n}']/(2*l(n))+0.5*[delp{n-1}*metricgradx{n}*delp{n-1}',delp{n-1}*metricgrady{n}*delp{n-1}',delp{n-1}*metricgradz{n}*delp{n-1}']/(2*l(n-1));   


% for i=1:n-1
%     delp{i}=[y(i+1,1)-y(i,1),y(i+1,2)-y(i,2),y(i+1,3)-y(i,3)];
% end
% delp{n}=[y(1,1)-y(n,1),y(1,2)-y(n,2),y(1,3)-y(n,3)];
% 
% for i=2:n-1
%     jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+0.5*[delp{i}*metricgradx{i}*delp{i}',delp{i}*metricgrady{i}*delp{i}',delp{i}*metricgradz{i}*delp{i}']/(2*l(i))+0.5*[delp{i-1}*metricgradx{i}*delp{i-1}',delp{i-1}*metricgrady{i}*delp{i-1}',delp{i-1}*metricgradz{i}*delp{i-1}']/(2*l(i-1));
% end
% jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+0.5*[delp{1}*metricgradx{1}*delp{1}',delp{1}*metricgrady{1}*delp{1}',delp{1}*metricgradz{1}*delp{1}']/(2*l(1))+0.5*[delp{n}*metricgradx{1}*delp{n}',delp{n}*metricgrady{1}*delp{n}',delp{n}*metricgradz{1}*delp{n}']/(2*l(n));   
% jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+0.5*[delp{n}*metricgradx{n}*delp{n}',delp{n}*metricgrady{n}*delp{n}',delp{n}*metricgradz{n}*delp{n}']/(2*l(n))+0.5*[delp{n-1}*metricgradx{n}*delp{n-1}',delp{n-1}*metricgrady{n}*delp{n-1}',delp{n-1}*metricgradz{n}*delp{n-1}']/(2*l(n-1));   


%% Jacobiandisp-jacobian for displacement produced by gait

parfor i=2:1:n-1
    jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),height(i,:),dimension);
end
jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),height(1,:),dimension);
jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),height(n,:),dimension);

%% Jacobianeqi-term that keeps points eqi distant from each other


parfor i=2:n-1;
    len=sqrt((y(i+1,:)-y(i-1,:))*((metric{i-1}+metric{i+1})/2)*(y(i+1,:)-y(i-1,:))');
    midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2))/2;
    betacos=(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i})/2))'/(l(i-1)*len);
    xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*l(i-1)*betacos/len;
    jacobianeqi(i,:)=midpoint-xhat;
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



% for i=2:n-1;
%     len=sqrt((y(i+1,:)-y(i-1,:))*(y(i+1,:)-y(i-1,:))');
%     l1=sqrt((y(i,:)-y(i-1,:))*(y(i,:)-y(i-1,:))');
%     midpoint=y(i-1,:)+(y(i+1,:)-y(i-1,:))/2;
%     betacos=(y(i+1,:)-y(i-1,:))*(y(i,:)-y(i-1,:))'/(l1*len);
%     xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*l1*betacos/len;
%     jacobianeqi(i,:)=midpoint-xhat;
% end
% 
%     len=sqrt((y(2,:)-y(n,:))*(y(2,:)-y(n,:))');
%     l1=sqrt((y(1,:)-y(n,:))*(y(1,:)-y(n,:))');
%     midpoint=y(n,:)+(y(2,:)-y(n,:))/2;
%     betacos=(y(2,:)-y(n,:))*(y(1,:)-y(n,:))'/(l1*len);
%     xhat=y(n,:)+(y(2,:)-y(n,:))*l1*betacos/len;
%     jacobianeqi(1,:)=midpoint-xhat;
% 
%     len=sqrt((y(1,:)-y(n-1,:))*(y(1,:)-y(n-1,:))');
%     l1=sqrt((y(2,:)-y(1,:))*(y(2,:)-y(1,:))');
%     midpoint=y(n-1,:)+(y(1,:)-y(n-1,:))/2;
%     betacos=(y(1,:)-y(n-1,:))*(y(n,:)-y(n-1,:))'/(l1*len);
%     xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*l1*betacos/len;
%     jacobianeqi(n,:)=midpoint-xhat;



%% Final gradient calculation
if invert == 0
    %totaljacobian=jacobiandisp/totalstroke-(((lineint)*jacobianstroke)/(1*totalstroke^2))/1+1*jacobianeqi;
    totaljacobian=jacobiandisp-(((lineint)*jacobianstroke)/(1*totalstroke))/1+1*jacobianeqi;
    dsdw=jacobiandisp;
    dy=dsdw(:);
else
    %totaljacobiantemp=jacobiandisp/totalstroke-(((lineint)*jacobianstroke)/(1*totalstroke^2))/1+1*jacobianeqi;
    totaljacobiantemp=jacobiandisp-(((lineint)*jacobianstroke)/(1*totalstroke))/1+1*jacobianeqi;
    for i=1:n
        totaljacobian(i,:)=totaljacobiantemp(n+1-i,:);
    end
end

%% minimizing negative of efficiency
f=-lineint/(totalstroke);
if nargout>1
    %g=-(jacobiandisp)-1*jacobianeqi;
    g=-totaljacobian(:);
end

%% Debugging and plotting

% for i=1:n
%     G(i)=y(i,1);
%     H(i)=y(i,2);
% %     P(i)=y(i,3);
%     I(i)=1*jacobianstroke(i,1);
%     J(i)=1*jacobianstroke(i,2);
% %     N(i)=1*jacobianstroke(i,3);
%     K(i)=jacobiandisp(i,1);
%     L(i)=jacobiandisp(i,2);
% %     Q(i)=jacobiandisp(i,3);
%     B(i)=totaljacobian(i,1);
%     C(i)=totaljacobian(i,2);
% %     D(i)=totaljacobian(i,3);
%     O(i)=jacobianeqi(i,1);
%     S(i)=jacobianeqi(i,2);
% %     R(i)=jacobianeqi(i,3);
% %     heightx(i)=height(i,1);
% %     heighty(i)=height(i,2);
% %     heightz(i)=height(i,3);
% %     heightcrosscheck(i)=height(i,:)*jacobiandisp(i,:)';
% end
% 
% 
% clf(figure(6)) %%jacobiandisp
% figure(6)
% scale=0;
% scale1=1;
% quiver(G,H,K,L,scale1)
% hold on
% % quiver(G,H,P,heightx,heighty,heightz,scale1)
% plot(G,H)
% axis equal
% hold off
% % 
% clf(figure(7)) %%jacobianstroke
% figure(7)
% scale=0;
% % quiver(G,H,K,L,scale)
% quiver(G,H,I,J,scale1)
% hold on
% plot(G,H)
% axis equal
% hold off
% % % % 
% clf(figure(8)) %%% totaljacobian
% figure(8)
% scale=0;
% quiver(G,H,B,C,scale1)
% hold on
% plot(G,H)
% axis equal
% hold off
% 
% clf(figure(3)) %%%jacobianeqi
% figure(3)
% scale1=1;
% quiver3(G,H,P,O,S,R,scale1)
% hold on
% plot3(G,H,P)
% axis equal
% hold off
% 
% 
% % figure(5)
% % scale=0;
% % quiver(G,H,B1,C1,scale1)
% % hold on
% % plot(G,H)
% % axis equal
% % hold off
% pause(0.1)
end

function a=jacobiandispcalculator3(p1,p2,p3,height,dimension)



% l=sqrt((p1(1)-p3(1))^2+(p1(2)-p3(2))^2+(p1(3)-p3(3))^2);
l1=0;
for i=1:1:dimension
    l1=l1+(p1(i)-p3(i))^2;
    base(1,i)=p3(i)-p1(i);
end
l=sqrt(l1);

for i=1:1:dimension
    jacobian(1,i)=0;
    perp1=zeros(1,dimension);
    perp1(i)=1;
    %parcomp=base*perp1'/norm(base);
    perp=perp1;%perp1-parcomp*base/norm(base);  %%recheck again
    for j=1:dimension-1
        for k=1:dimension-j
            veca=zeros(1,dimension);
            vecb=zeros(1,dimension);
            veca(j)=1;
            vecb(j+k)=1;
            f=(j-1)*dimension-(j*(j-1))/2+k;
            jacobian(1,i)=jacobian(1,i)+0.5*height(f)*((veca*perp')*(vecb*base')-(vecb*perp')*(veca*base'));
        end
    end
end

a=jacobian;

% %% local coordinates
% xparallel1=[p3(1)-p1(1),p3(2)-p1(2),p3(3)-p1(3)];
% xparallel=xparallel1/norm(xparallel1);
% 
% xneeded1=[p2(1)-p1(1),p2(2)-p1(2),p2(3)-p1(3)];
% xneeded=xneeded1/norm(xneeded1);
% 
% xper21=cross(xparallel,xneeded);
% xper2=xper21/norm(xper21);
% 
% xper11=cross(xper2,xparallel);
% xper1=xper11/norm(xper11);
% 
% %% jacobian calculation
% per1flux=-xper2;
% per2flux=xper1;
% 
% heightper1=height*per1flux';
% heightper2=height*per2flux';
% 
% jacobian=0.5*l*heightper1*xper1+0.5*l*heightper2*xper2;
% a=jacobian;

end

function [net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,tspan,ConnectionEval,IntegrationMethod,resolution)
% Evaluate the displacement and cost for the gait specified in the
% structure GAIT when executed by the system described in the structure
% S.
%
% S should be a sysplotter 's' structure loaded from a file
% sysf_FILENAME_calc.mat (note the _calc suffix)
%
% P should be a structure with fields "phi" and "dphi", returning a
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
            startshape = p.phi(0);
            startshapelist = num2cell(startshape);
            beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'cubic');
            net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
                -sin(beta_theta) cos(beta_theta) 0;...
                0 0 1]*net_disp_orig;
            
        otherwise
			error('Unknown method for integrating motion');
	end

	
	% Convert the final motion into its representation in optimal
	% coordinates
	startshape = p.phi(0);
	startshapelist = num2cell(startshape);
	beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'cubic');
	net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
		-sin(beta_theta) cos(beta_theta) 0;...
		0 0 1]*net_disp_orig;

	
end

% Evaluate the body velocity and cost velocity (according to system metric)
% at a given time
function [xi, dcost] = get_velocities(t,s,gait,ConnectionEval)

	% Get the shape and shape derivative at the current time
	shape = gait.phi(t);
	shapelist = num2cell(shape);
	dshape = gait.dphi(t);
	
	
	% Get the local connection and metric at the current time, in the new coordinates	
	switch ConnectionEval
		case 'functional'
			
			A = s.A_num(shapelist{:})./s.A_den(shapelist{:});
			
			M = s.metric(shapelist{:});

		case 'interpolated'
			
			A = -cellfun(@(C) interpn(s.grid.eval{:},C,shapelist{:}),s.vecfield.eval.content.Avec);
			
			M =  cellfun(@(C) interpn(s.grid.metric_eval{:},C,shapelist{:}),s.metricfield.metric_eval.content.metric);
			
		otherwise
			error('Unknown method for evaluating local connection');
	end
	
	% Get the body velocity at the current time
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
