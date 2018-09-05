function [X, Y] = squashed_rectangle(w,h,cap_fraction,points)
% Get the X,Y points on a "squashed rectangle" as defined in the OmniGraffle
% drawing program (a rectangle for which cap_fraction of the width is a split
% ellipse. To visually match OmniGraffle figures, cap_fraction should be .2
%
% All the points returned will be on the elliptical end caps, with the
% final points on the caps the endpoints of the straight edges. the points
% parameter specifies how many points are included in each half-ellipse.

%Create a parameterization for the points on the ellipse
t = linspace(0,pi,points)';

%Generate half of 'an ellipse with height h and width cap_percent*w'
a = h/2;
b = cap_fraction*w/2;

X = -b*sin(t);
Y = a*cos(t);

%Offset the cap
X = X-(1-cap_fraction)*w/2;

%Mirror the cap
X = [X;-X;X(1)];
Y = [Y;-Y;Y(1)];

end