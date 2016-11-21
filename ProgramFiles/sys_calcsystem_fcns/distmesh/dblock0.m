function d=dblock0(p,x1,x2,y1,y2,z1,z2)

% Distances to each face

df1=x1-p(:,1); % X distance from -x face
df2=-x2+p(:,1); % X distance from +x face 
df3=y1-p(:,2);  % Y distance from -y edge
df4=-y2+p(:,2); % Y distance from +y edge
df5=z1-p(:,3); % Z distance from -z face
df6=-z2+p(:,3); % Z distance from +z face 


% Distances to each edge

d5=sqrt(d1.^2+d3.^2); % L2 distance from -x-y corner
d6=sqrt(d1.^2+d4.^2); % L2 distance from +x-y corner
d7=sqrt(d2.^2+d3.^2); % L2 distance from -x+y corner
d8=sqrt(d2.^2+d4.^2); % L2 distance from +x+y corner

d=-min(min(min(-d1,-d2),-d3),-d4);  
% (distance from closest y edge)
% (distance from closest edge)

% Correct for other regions

ix=d1>0 & d3>0;  % if in -- corner, use squared distance
d(ix)=d5(ix);
ix=d1>0 & d4>0;
d(ix)=d6(ix);
ix=d2>0 & d3>0;
d(ix)=d7(ix);
ix=d2>0 & d4>0;
d(ix)=d8(ix);