function s = softspace(startPoint, endPoint, varargin)
%more linear version of cosspace

%set the default spacing if not provided
if ~isempty(varargin)
	n_points = varargin{1};
else
	n_points = 100;
end


%Build the acceleration and decelleration curves
c = cosspace(0,1,floor(n_points/2)-mod(floor(n_points/2),2));

%Build the linear section
l = linspace(0,pi/2,ceil(n_points/2)+1+mod(floor(n_points/2),2));

%Merge the two sections together
s = [c(1:length(c)/2) l(2:end-1)+c(length(c)/2) c(length(c)/2:end)+l(end)];

%Scale the output
s = (s/(1+pi/2) * (endPoint-startPoint)) + startPoint;