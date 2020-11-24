function [zcg,t,n,del] = ellipse(a,b,npts)
% -----------------
% E Kanso, 14 april 2004


% -----------------INPUT
% 
% a & b   major and minor axes of the ellipse
%
% npts    number of panels used to discretize the ellipse
%
% ----------------

% -----------------OUTPUT
% 
% zcg    position of collocation pts w.r.t inertial 
%        frame attached at the c.o.m of the ellipse
%
% t      components of vectors tangent to panels
%
% n      components of outward normal vectors
%
% del    panel length
%
% ----------------

beta  = linspace(0,2*pi,npts+1);
alpha =  beta';

Px = a*cos(alpha);  Py = b*sin(alpha);
Px = Px(npts+1:-1:1); Py = Py(npts+1:-1:1);

% control pts
xcg = (Px(1:1:npts) + Px(2:1:npts+1))/2; 
ycg = (Py(1:1:npts) + Py(2:1:npts+1))/2;
zcg = [xcg,ycg];

% Panel length
Xrel =  diff(Px); Yrel =  diff(Py);   
del = sqrt(Xrel.^2+Yrel.^2);

% Tangent and normal vectors to panels
tx = Xrel./del; ty = Yrel./del;
t = [tx,ty]; n = [-ty,tx];

