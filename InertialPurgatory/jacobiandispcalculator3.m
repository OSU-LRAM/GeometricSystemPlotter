function a=jacobiandispcalculator3(p1,p2,p3,height)

l = sqrt((p1(1)-p3(1))^2+(p1(2)-p3(2))^2);%+(p1(3)-p3(3))^2);

%% local coordinates
xparallel1 = [p3(1)-p1(1),p3(2)-p1(2)];%,p3(3)-p1(3)];
xparallel = [xparallel1/norm(xparallel1) 0];

xneeded1 = [p2(1)-p1(1),p2(2)-p1(2)];% ,p2(3)-p1(3)];
xneeded = [xneeded1/norm(xneeded1) 0];

xper21 = cross(xparallel,xneeded);
xper2 = [xper21/norm(xper21)];

xper11 = cross(xper2,xparallel);
xper1 = [xper11/norm(xper11)];

%% jacobian calculation
per1flux = -xper2;
per2flux = xper1;

heightper1 = height*per1flux';
% heightper2=height*per2flux';

jacobian = 0.5*l*heightper1*xper1; %+0.5*l*heightper2*xper2;
jacobian(3) = [];
a = jacobian;


% l1=0;
% for i=1:1:dimension
%     l1=l1+(p1(i)-p3(i))^2;
%     base(1,i)=p3(i)-p1(i);
% end
% l=sqrt(l1);
% 
% for i=1:1:dimension
%     jacobian(1,i)=0;
%     perp1=zeros(1,dimension);
%     perp1(i)=1;
%     %parcomp=base*perp1'/norm(base);
%     perp=perp1;%perp1-parcomp*base/norm(base);  %%recheck again
%     for j=1:dimension-1
%         for k=1:dimension-j
%             veca=zeros(1,dimension);
%             vecb=zeros(1,dimension);
%             veca(j)=1;
%             vecb(j+k)=1;
%             f=(j-1)*dimension-(j*(j-1))/2+k;
%             jacobian(1,i)=jacobian(1,i)+0.5*height(f)*((veca*perp')*(vecb*base')-(vecb*perp')*(veca*base'));
%         end
%     end
% end
% 
% a=jacobian;

end
