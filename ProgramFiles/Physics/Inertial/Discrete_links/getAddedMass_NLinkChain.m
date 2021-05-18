function [FullMassMatrix,J,J_Full,h,local_inertias] = getAddedMass_NLinkChain(shape,s)
% Nathan Justus, 5/25/2020
% Adaptation of work by Eva Kanso 2004 for sysplotter

%-----------INPUTS------------
% shape - column vector of shape variables
% s - structure containing geometry and physical properties of N-link chain

%----------OUTPUT-------------
% AddedMassMatrix - locked inertia matrix (size: 3x3)


%Normalize swimmer to length 1

%Number of elliptic links in chain geometry
numLinks = numel(s.geometry.linklengths);
if isfield(s.geometry,'linkSeparation')
    linkSep = s.geometry.linkSeparation;
else
    linkSep = 0;
    s.geometry.linkSeparation = 0;
end
totalLength = sum(s.geometry.linklengths)+linkSep*(numLinks-1);
s.geometry.length = s.geometry.length/totalLength;
s.geometry.linklengths = s.geometry.linklengths/totalLength;
s.geometry.linkSeparation = s.geometry.linkSeparation/totalLength;

%If not specified, turn 'on' all links
if ~isfield(s.geometry,'activelinks')
    s.geometry.activelinks = ones(1,numLinks);
end

%If not specified, turn off inter-link interaction to save computation time
if ~isfield(s.physics,'interaction')
    s.physics.interaction = 'off';
end


%Build geometry of body shape
[zc,t,n,del,orientations,linkCenters,h] = buildNLinkBody(shape,s);

%Get full mass matrix and jacobians
[FullMassMatrix,J,J_Full,local_inertias] = fullMassMatrix_NLinkChain_byLink(zc,t,n,del,s,shape,orientations,linkCenters);

end

function [allPanelPoints,t,n,del,orientations,zg,h] = buildNLinkBody(shape,s)

% Nathan Justus 5/26/20

% Generates panel positions, normal/tangent vectors, and other relevant 
% spatial information based on given geometry

%-----------INPUTS-------------
% shape - column vector of shape variables
% s - structure containing geometry and physical properties of N-link chain

%-----------OUTPUTS------------
% allPanelPoints - position of ellipse points w.r.t base frame from s
% t - components of normal vectors tangent to ellipse panels
% n - components of outward normal vectors for ellipse panels
% del - ellipse panel lengths
% orientations - vector of ellipse orientations w.r.t. swimmer frame
% zg - matrix of ellipse cg's w.r.t. swimmer frame
% h - structure encoding swimmer configuration

%Deal joint angles depending on if we're using mode shapes
if isfield(s.geometry,'modes')
    jointAngles = s.geometry.modes*shape;
else
    jointAngles = shape;
end

% Number of elliptic links in chain geometry
numLinks = numel(s.geometry.linklengths);
% List of ellipse lengths
linklengths = s.geometry.linklengths;
% Define ratio parameters for ellipse:
aspectRatio = s.geometry.link_shape_parameters{1}.aspect_ratio;
% Number of panels comprising ellipse
npts = 100;  
% Separation distance between ellipses at joints - avoids singularities
linkSep = s.geometry.linkSeparation/2;

% Check to make sure that the given information makes sense with geometry
if numel(jointAngles)+1 ~= numLinks
    error('Wrong number of joints for this system');
end

%XY locations of c.o.m. for each ellipse w.r.t. first link c.o.m.
zg = zeros(numLinks,2);
%Orientation of each link
orientation = 0;
%Get cg of each link w.r.t. first link cg
for i = [1:numLinks-1]
    %Get semi-major axis of current ellipse and next ellipse
    a1 = linklengths(i)/2;
    a2 = linklengths(i+1)/2;
    %Grab cg of current ellipse
    lastCG = zg(i,:);
    %Move to point at right edge of this ellipse (joint location)
    rightEdge = lastCG + [(a1+linkSep)*cos(orientation),(a1+linkSep)*sin(orientation)];
    %If there are an even number of links and this joint is the middle,
    if mod(numLinks,2) == 0 && i == numLinks/2
        %Log middle mosition as joint location
        middlePosition = rightEdge;
    end
    %Rotate by the joint angle
    orientation = orientation + jointAngles(i);
    %Move to cg of new link
    zg(i+1,:) = rightEdge + [(a2+linkSep)*cos(orientation),(a2+linkSep)*sin(orientation)];
end

%Deal out link orientations with respect to middle link
orientations = [];
%If odd number of links, frame alighned with middle link frame
if mod(numLinks,2) == 1
    for i = [1:numLinks]
        if i < numLinks/2
            orientations(i) = -sum(jointAngles(i:floor(numLinks/2)));
        elseif i == floor(numLinks/2)+1
            orientations(i) = 0;
        elseif i > numLinks/2 + 1
            orientations(i) = sum(jointAngles(floor(numLinks/2)+1:i-1));
        end
    end
%If even number of links, frame is oriented at average of the two
%center-most link orientations
else
    orientations(1) = 0;
    for i = [1:numel(jointAngles)]
        orientations(i+1) = orientations(i) + jointAngles(i);
    end
    %Orientation of base frame with respect to leftmost link
    middleOrientation = mean(orientations(numLinks/2:numLinks/2+1));
    orientations = orientations - middleOrientation;
end

%Move/rotate center of masses for each link based off base frame
if s.geometry.baseframe == 'center'
    %If odd number of links, get orientation and frame location from center
    %ellipse
    if mod(numLinks,2) == 1
        %Get orientation of middle link w.r.t. left link
        middleOrientation = sum(jointAngles(1:numel(jointAngles)/2));
        %Get cg position of middle link w.r.t. left link cg
        middlePosition = zg(floor(numLinks/2)+1,:);
    end
    %Shift cgs so frame center is at 0,0
    zg = zg - repmat(middlePosition,numLinks,1);
    %Rotate cgs about center link so frame is aligned with working x,y
    R = [cos(-middleOrientation),-sin(-middleOrientation);...
        sin(-middleOrientation),cos(-middleOrientation)];
    zg = (R*zg')';
end

%Make structure with link cg positions and lengths for logging
h.pos = [zg,orientations'];
h.lengths = s.geometry.linklengths';

%Move ellipses to correct cgs/orientations
allPanelPoints = [];
t = [];
n = [];
del = [];
%For each link in the ellipse
for i = [1:numLinks]
    %If that link is turned 'on'
    if s.geometry.activelinks(i)
        %Get ellipse data
        a = linklengths(i)/2;
        b = a*aspectRatio;
        [zcg,ts,ns,dels] = ellipse(a,b,npts);
        %Make rotation matrix
        R = [cos(orientations(i)),-sin(orientations(i));...
            sin(orientations(i)),cos(orientations(i))];
        %Rotate and translate ellipse to correct position/orientation
        thisPanel = (R*zcg')' + zg(i,:);
        allPanelPoints = [allPanelPoints;thisPanel];
        %Rotate tangent vectors
        thisT = (R*ts')';
        t = [t;thisT];
        %Get normal vectors
        thisN = (R*ns')';
        n = [n;thisN];
        %Get panel lengths
        del = [del;dels];
    end
end


end

function [An,Bt,Phi] = influence_NLinkChain_byLink(zc,t,n,del,N)
% -----------------
% E Kanso, 14 april 2004

% -----------------INPUT
% 
% zc     position of collocation pts 
% t      components of vectors tangent to panels
% n      components of outward normal vectors 
% del    panel length
% N      number of panels
%
% zc, t and n  are w.r.t base frame
%
% ----------------


% -----------------INTERNAL VARIABLES
% 
% Xrel(i,j) & Yrel(i,j)    coordinates of control pt i (Ci) relative 
%                          to control pt j (Cj) w.r.t inertial frame
%
% Cn(i,j) & Ct(i,j)        normal and tangential coordinates of Ci relative
%                          to Cj w.r.t. a frame attached to the panel j
%                           
% Vn(i,j) & Vt(i,j)         normal and tangential velocities induced 
%                           at Ci due to a constant source distribution 
%                           at panel j,  w.r.t. a frame attached to panel j
%                           
% Vx(i,j) & Vy(i,j)         velocities Vn(i,j) and Vt(i,j)
%                           expressed w.r.t. inertial frame
%
% ----------------


% -----------------OUTPUT 
% 
% An(i,j) & Bt(i,j)       normal and tangential velocities induced at Ci
%                         due to a constant source distribution at panel j 
%                         Expressed w.r.t a frame attached to the panel i
%
% Phi                     potential function
%
% ----------------


% initialize
Xc = zeros(N,N); Yc = zeros(N,N);
tx = zeros(N,N); ty = zeros(N,N);
nx = zeros(N,N); ny = zeros(N,N);
Ds = zeros(N,N);


% assign
Xc = zc(:,1)*ones(1,N); Yc = zc(:,2)*ones(1,N);
tx = t(:,1)*ones(1,N);  ty = t(:,2)*ones(1,N);
nx = n(:,1)*ones(1,N);  ny = n(:,2)*ones(1,N);

Ds = 0.5.*ones(N,1)*del';


% compute Xrel(i,j) and Yrel(i,j) 
XREL = Xc-Xc';  
YREL = Yc-Yc';  

% compute Cn(i,j) and Ct(i,j)
Cn = XREL.*nx' + YREL.*ny';
Ct = XREL.*tx' + YREL.*ty';


% compute Vn(i,j) and Vt(i,j) 
temp1 = Cn.^2;
temp2 = (Ct + Ds).^2;
temp3 = (Ct - Ds).^2;
Vt_Num = temp2 + temp1;
Vt_Den = temp3 + temp1;
Vt = log(Vt_Num./Vt_Den);
Vn_Num = 2.*Cn.*Ds;
Vn_Den = Ct.^2 + temp1 - Ds.^2;
angle = 2.*atan(Vn_Num./Vn_Den); 
Vn = angle + 2*pi*eye(N,N);

% compute Vx(i,j) and Vy(i,j) 
Vx =   Vn.*nx' + Vt.*tx';
Vy =   Vn.*ny' + Vt.*ty';
        
% compute An(i,j) and Bt(i,j) 
An =  Vx.*nx + Vy.*ny;
Bt =  Vx.*tx + Vy.*ty;
        
% compute Phi
Phi = - Ct.*Vt - Cn.*Vn - Ds.*log((temp2 + temp1).*(temp3 + temp1));

end

function linkAddedMass = getAddedMassByLink(s,zcg,t,n,del,linkNumber,inv_An,Phi)

% Nathan Justus 5/26/20

% Generates added mass for a given link in swimmer frame

%-----------INPUTS-------------
% s - structure containing geometry and physical properties of N-link chain
% zcg - N by 2 matrix of panel locations where N is # of panels
% t - N by 2 matrix of panel unit tangent vectors
% n - N by 2 matrix of panel unit normal vectors
% del - N-length vector of panel lengths
% linkNumber - scalar identifier for ellipse, relative to ~all~ links, l->r
% inv_An - inverse of induced normal velocities - see influence function
% Phi - potential function - see influence function for details

%-----------OUTPUTS------------
% linkAddedMass - 3 by 3 added mass matrix for chosen link in swimmer frame

%Get system information
%Density of fluid swimmer is submerged in
density = s.physics.fluid_density;
%Number of links turned 'on' in system
numLinks = sum(s.geometry.activelinks);
%Redo linkNumber to be relative to only links turned 'on'
linkNumber = sum(s.geometry.activelinks(1:linkNumber));
%number of panels in body
totalPoints = numel(del);
%number of panels per ellipse
npts = totalPoints/numLinks;

%normal panel velocity due to angular rotation at unit speed
avel = zcg(:,1).*n(:,2) - zcg(:,2).*n(:,1);
%normal panel velocity due to unit velocity in x direction (right)
xvel = ones(totalPoints,1).*n(:,1);
%normal panel velocity due to unit velocity in y direction (up)
yvel = ones(totalPoints,1).*n(:,2);

%For every active link
for i = [1:numLinks]
    
    %Body velocity components - zero everywhere but at corresponding links
    vfn_x = zeros(totalPoints,1);
    vfn_y = zeros(totalPoints,1);
    vfn_a = zeros(totalPoints,1);
    vfn_x(npts*(i-1)+1:npts*i) = xvel(npts*(i-1)+1:npts*i);
    vfn_y(npts*(i-1)+1:npts*i) = yvel(npts*(i-1)+1:npts*i);
    vfn_a(npts*(i-1)+1:npts*i) = avel(npts*(i-1)+1:npts*i);
    vfn{i}.x = vfn_x;
    vfn{i}.y = vfn_y;
    vfn{i}.a = vfn_a;
    
    %Source density distributions
    sigma{i}.x = inv_An*vfn_x;
    sigma{i}.y = inv_An*vfn_y;
    sigma{i}.a = inv_An*vfn_a;
    
    %Potential functions
    phi{i}.x = Phi*sigma{i}.x;
    phi{i}.y = Phi*sigma{i}.y;
    phi{i}.a = Phi*sigma{i}.a;

end

%Storage for added mass contributions
mass = zeros(3);
%Compute nondimensionalized added masses for each active ellipse's effect 
%on every other active ellipse
i = linkNumber;
for j = [1:numLinks]
    %Added mass of ellipse just by itself
    if i == j
        xx = density*sum(phi{i}.x.*vfn{i}.x.*del);
        yy = density*sum(phi{i}.y.*vfn{i}.y.*del);
        aa = density*sum(phi{i}.a.*vfn{i}.a.*del);
        xy = density*sum(phi{i}.x.*vfn{i}.y.*del);
        xa = density*sum(phi{i}.x.*vfn{i}.a.*del);
        ya = density*sum(phi{i}.y.*vfn{i}.a.*del);
        thisMass = [xx,xy,xa;...
                    xy,yy,ya;...
                    xa,ya,aa];
        %Add contribution to total added mass matrix
        mass = mass+thisMass;
    %Added mass of ellipse from interaction with other ellipse
    else
        xx = density*sum(phi{i}.x.*vfn{j}.x.*del);
        yy = density*sum(phi{i}.y.*vfn{j}.y.*del);
        aa = density*sum(phi{i}.a.*vfn{j}.a.*del);
        xy = .5*density*(sum(phi{i}.x.*vfn{j}.y.*del)+sum(phi{j}.x.*vfn{i}.y.*del));
        xa = .5*density*(sum(phi{i}.x.*vfn{j}.a.*del)+sum(phi{j}.x.*vfn{i}.a.*del));
        ya = .5*density*(sum(phi{i}.y.*vfn{j}.a.*del)+sum(phi{j}.y.*vfn{i}.a.*del));
        thisMass = [xx,xy,xa;...
                    xy,yy,ya;...
                    xa,ya,aa];
        %Add contribution to total added mass matrix
        mass = mass+thisMass;
    end
end

%Return total added mass matrix
linkAddedMass = mass;

end
      
function [FullMassMatrix,jacobians,jacobians_full,local_inertias] = fullMassMatrix_NLinkChain_byLink(zcg,t,n,del,s,shape,orientations,linkCenters)
% Nathan Justus 5/26/20

% Generates full mass matrix for system in swimmer frame

%-----------INPUTS-------------
% zcg - N by 2 matrix of panel locations where N is # of panels
% t - N by 2 matrix of panel unit tangent vectors
% n - N by 2 matrix of panel unit normal vectors
% del - N-length vector of panel lengths
% s - structure containing geometry and physical properties of N-link chain
% shape - shape variables for swimmer
% orientations - vector of ellipse orientations w.r.t. swimmer frame
% linkCenters - matrix of ellipse cg's w.r.t. swimmer frame

%-----------OUTPUTS------------
% FullMassMatrix - (3+S) by (3+S) mass matrix where S is # of shape vars
% jacobians - M-length cell vector of 3 by S matrices of link body
%             jacobians where M is number of system ellipses
% jacobians_full - M-length cell vector of 3 by (3+S) matrices 
%                  of link+spatial body jacobians

%Storage for mass matrix
FullMassMatrix = zeros(3+numel(shape));
%Number of links in system
numLinks = size(linkCenters,1);
%Density of surrounding fluid
density = s.physics.fluid_density;
%Number of panels in system
totalPoints = numel(del);
%Initialize cell storage for jacobians
jacobians = {};
jacobians_full = {};
%Initialize cell storage for local inertias
local_inertias = {};

%If inter-panel interaction is turned on
if strcmp(s.physics.interaction,'on')
    %Get panel relative influences and potential function
    [An,Bt,Phi] = influence_NLinkChain_byLink(zcg,t,n,del,totalPoints);
    %Invert induced panel normal velocities
    inv_An = inv(An);
end

%Aspect ratio of ellipses in chain
aspectRatio = s.geometry.link_shape_parameters{1}.aspect_ratio;

%For each ellipse in the chain
for link = [1:numLinks]
    
    %Initialize local inertia as zero mass matrix
    local_inertias{link} = zeros(3);
    
    %Compute jacobian for corresponding link in chain
    thisJacobian = getJacobian_NLinkChain(s,shape,link,orientations,linkCenters);
    jacobians_full{link} = thisJacobian;
    jacobians{link} = thisJacobian(:,4:end);
    
    %If the link is turned 'on'
    if s.geometry.activelinks(link)
        
        %Get ellipse parameters for this link
        linkLength = s.geometry.linklengths(link);
        a = linkLength/2;
        b = a*aspectRatio;
        
        %Calculate ellipse mass in link frame
        mass_link = pi*a*b;
        rotational_inertia_link = mass_link*(a^2+b^2)/4;
        standardmass = diag([mass_link,mass_link,rotational_inertia_link]);

        %If interaction is turned on, use potential field to get added mass
        if strcmp(s.physics.interaction,'on')
            %Get link added mass in base swimmer frame
            linkAddedMass = getAddedMassByLink(s,zcg,t,n,del,link,inv_An,Phi);
            %Push forward calculated mass from base frame to link frame
            pushforwardJ = thisJacobian(:,1:3);
            iJ = inv(pushforwardJ);
            iJt = inv(pushforwardJ');
            linkAddedMass = iJt*linkAddedMass*iJ;
        %If interaction is turned off, just use known ellipse added mass
        else
            m_xx = density*pi*b^2;
            m_yy = density*pi*a^2;
            m_aa = density*pi/8*(a^2-b^2)^2;
            linkAddedMass = diag([m_xx,m_yy,m_aa]);
        end

        %Log local inertia in link frame
        local_inertias{link} = linkAddedMass + standardmass;
        
        %Pull back added mass + standard mass to base frame through full joint jacobians
        thisFullMass = thisJacobian'*(linkAddedMass+standardmass)*thisJacobian;
        %Add this to the total mass matrix
        FullMassMatrix = FullMassMatrix+thisFullMass;
    end

end

end

function bodyJacobian = getJacobian_NLinkChain(s,shape,linkNumber,orientations,linkCenters)
% Nathan Justus 5/26/20

% Generates body jacobian for link wrt swimmer frame

%-----------INPUTS-------------
% s - structure containing geometry and physical properties of N-link chain
% shape - shape variables for swimmer
% linkNumber - identifier for link in chain wrt ~all~ other links
% orientations - vector of ellipse orientations w.r.t. swimmer frame
% linkCenters - matrix of ellipse cg's w.r.t. swimmer frame

%-----------OUTPUTS------------
% bodyJadobian - 3 by (3+S) body jacobian matrix for given link wrt
%                swimmer base frame where S is # of shape variables

%If using shape modes
if isfield(s.geometry,'modes')
    %Import them to convert shape variables to joint angles
    jointAngles = s.geometry.modes*shape;
%Otherwise
else
    %Just use joint angles directly
    jointAngles = shape;
end

%Get distance between edge of ellipse and actual joint location
linkSeparation = s.geometry.linkSeparation/2;

%If the base frame is in the middle of the center link
if s.geometry.baseframe == 'center'
    %Figure out which link number corresponds to center link
    midLink = numel(jointAngles)/2+1;
    %Initialize local body jacobian
    bodyJacobian = zeros(3,numel(jointAngles));
    %Initialize adjoint chain
    adjoints = eye(3);
    %If the chosen link is to the right of the center link, we will walk
    %left along the chain back to the center
    if linkNumber > midLink
        for i = [linkNumber:-1:midLink+.5]
            linkWidth1 = s.geometry.linklengths(i)/2;
            linkStep1 = linkWidth1+linkSeparation;
            linkWidth2 = s.geometry.linklengths(i-1)/2;
            linkStep2 = linkWidth2+linkSeparation;
            %Step halfway along a link
            stepRight1 = invAdj2D([linkStep1,0,0]);
            stepRight2 = invAdj2D([linkStep2,0,0]);
            %If there are an even number of ellipses and almost at middle
            if floor(midLink) ~= midLink && i-midLink < 1
                %Rotate frame by half joint angle
                rotate = invAdj2D([0,0,jointAngles(i-1)/2]);
                %Get the jacobian column for a half-angle at joint location
                jacobianCol = adjoints*stepRight1*rotate*[0;0;.5];
                bodyJacobian(:,i-1) = jacobianCol;
                %Move adjoint chain to base frame
                adjoints = adjoints*stepRight1*rotate;
            %If not to middle yet or for an odd number of ellipses
            else
                %Rotate frame by the joint angle
                rotate = invAdj2D([0,0,jointAngles(i-1)]);
                %Get the jacobian column from chain so far
                jacobianCol = adjoints*stepRight1*rotate*[0;0;1];
                bodyJacobian(:,i-1) = jacobianCol;
                %Move adjoint chain to center of proximal link for next step
                adjoints = adjoints*stepRight1*rotate*stepRight2;
            end
        end
    %If the chosen link is to the left of the center link, we will walk
    %right along the chain back to the center
    elseif linkNumber < midLink
        for i = [linkNumber:midLink-.5]
            linkWidth1 = s.geometry.linklengths(i)/2;
            linkStep1 = linkWidth1+linkSeparation;
            linkWidth2 = s.geometry.linklengths(i+1)/2;
            linkStep2 = linkWidth2+linkSeparation;
            %Step halfway along a link to the left
            stepLeft1 = invAdj2D([-linkStep1,0,0]);
            stepLeft2 = invAdj2D([-linkStep2,0,0]);
            %If even number of links and almost to middle
            if floor(midLink) ~= midLink && midLink-i < 1
                %Rotate frame opposite of half joint angle, since moving backwards
                rotate = invAdj2D([0,0,-jointAngles(i)/2]);
                %Get the jacobian column for a half-angle at joint location
                jacobianCol = adjoints*stepLeft1*rotate*[0;0;-.5];
                bodyJacobian(:,i) = jacobianCol;
                %Move adjoint chain to center of base frame
                adjoints = adjoints*stepLeft1*rotate;
            else
                %Rotate frame opposite of joint angle, since moving backwards
                rotate = invAdj2D([0,0,-jointAngles(i)]);
                %Get jacobian column from chain so far
                jacobianCol = adjoints*stepLeft1*rotate*[0;0;-1];
                bodyJacobian(:,i) = jacobianCol;
                %Move adjoint chain to center of proximal link for next step
                adjoints = adjoints*stepLeft1*rotate*stepLeft2;
            end
        end
    end

    %If using shape modes, pull jacobian through them to get jacobian wrt
    %actual shape variables instead of joint angles
    if isfield(s.geometry,'modes')
        bodyJacobian = bodyJacobian*s.geometry.modes;
    end
    
    %Tack on 3x3 adjoints matrix to account for body movements
    bodyJacobian = [adjoints,bodyJacobian];  
end

end        
        
function iadj = invAdj2D(g)
% Nathan Justus 5/26/20

% Get inverse adjoint representation of a rotation/translation

%-----------INPUTS-------------
% g - list holding x,y,alpha specifying translation/rotation amount

%-----------OUTPUTS------------
% iadj - inverse adjoint operator for body jacobian construction

%Deal out variables
x = g(1);
y = g(2);
a = g(3);

%Construct adjoint operator
adj = [cos(a),-sin(a),y;sin(a),cos(a),-x;0,0,1];
%Take the inverse
iadj = inv(adj);

end   