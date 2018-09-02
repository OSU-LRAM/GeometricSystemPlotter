function [B,h,J, J_full] = fat_chain_from_links(geometry, display, jointangles)
% Generate a set of rounded-end links based on link lengths and joint angles in
% a kinematic chain
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are
%           specified in the documentation of N_link_chain.
%
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   display: Structure containing information about how the system should
%       be illustrated. Field required by this function is:
%
%       aspect: ratio of width to length of the system
%
%
%   jointangles: A vector of the angles between the links. Must be one
%       element shorter than the linklengths vector


    % Get the positions and scaled lengths of the links
    [h,J, J_full] = N_link_chain(geometry,jointangles);
    
    % Convert the backbone locations into SE(2) transforms
    g = vec_to_mat_SE2(h.pos);
    
    %%%%%%%%
    % Create the links
    
    points = 50;             % Number of points per link
    N_links = size(g,3);      % Number of links in the chain
    B = [];                   % Array to hold points on chain
    
    % Iterate over links
    for idx = 1:N_links
        
        %%%
        % Create an ellipse-ended rectangle for the link 
        
        % Make the shape slightly shorter than this link (so that there is
        % a whitespace gap between them when plotted)
        length = h.lengths(idx)/1.05; %
        
        % Generate the width of the shape from the aspect ratio multiplied
        % by the total length of the system
        width = display.aspect_ratio * sum(h.lengths);
        
        % get the points on the ellipse-ended rectangle
        C = squashed_rectangle(length,width,0.2,points);
        
        % Transform rectangle to link position
        C = g(:,:,idx) * C;
        
        % Append C onto the end of B, and add a NaN so that the links will
        % plot discontinuously
        B = [B, C, NaN(3,1)];
        
    end

end


function XY = squashed_rectangle(length,width,cap_fraction,points)

    % Get the X,Y points on a "squashed rectangle" as defined in the OmniGraffle
    % drawing program (a rectangle for which cap_fraction of the width is a split
    % ellipse. To visually match OmniGraffle figures, cap_fraction should be .2
    %
    % All the points returned will be on the elliptical end caps, with the
    % final points on the caps the endpoints of the straight edges. the points
    % parameter specifies how many points are included in each half-ellipse.

    %Create a parameterization for the points on the ellipse
    t = linspace(0,pi,points)';

    %Generate an ellipse whose 'a' axis is the width of the link, and whose
    % 'b' axis is the cap fraction times the link length 
    a = width/2;
    b = cap_fraction*length/2;

    X = -b*sin(t);
    Y = a*cos(t);

    %Offset the cap
    X = X-(1-cap_fraction)*length/2;

    %Mirror the cap
    X = [X;-X;X(1)];
    Y = [Y;-Y;Y(1)];
    
    % Augment these values with a column of ones to make them amenable to
    % SE(2) transformations
    XY = [X Y ones(2*points+1,1)]';

end