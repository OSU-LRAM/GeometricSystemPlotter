function a=jacobiandispcalculator3(p1,p2,p3,ccf,dimension)
%%%%%%%%%
%
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. 
% Its input arguments are the coordinates of the (i-1)th, ith and (i+1)th point,     
% CCF value at point i(ccf) and the dimension of the shape space (dimension)
%
%%%%%%%%%

l1=0; % variable for calculating length of the line segment joining the (i-1)th point with the (i+1)th point
base = zeros(1,dimension);
for i=1:numel(base)
    l1=l1+(p1(i)-p3(i))^2;
    base(1,i)=p3(i)-p1(i); % vector connecting the (i-1)th point and (i+1)th point  
end
%l=sqrt(l1); % length of the line segment joining the (i-1)th point with the (i+1)th point

jacobian = zeros(1,dimension);
for i=1:dimension
%    jacobian(1,i)=0;
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