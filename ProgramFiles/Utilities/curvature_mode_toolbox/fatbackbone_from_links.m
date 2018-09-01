function [B,h,J, J_full] = fatbackbone_from_links(linklengths,jointangles,L,baseframe,width)
% Draw a set of squashed-rectangle links

    % Get the positions and scaled lengths of the links
    [h,J, J_full] = backbone_from_links(linklengths,jointangles,L,baseframe);
    
    % Convert the backbone locations into SE(2) transforms
    g = vec_to_mat_SE2(h.pos);
    
    %%%%%%%%
    % Create the links
    
    points = 50;             % Number of points per link
    N_links = size(g,3);      % Number of links in the chain
    B = [];                   % Array to hold points on chain
    
    % Iterate over links
    for idx = 1:N_links
        
        % Create a squashed rectangle slightly shorter than this link (so
        % that there is a whitespace gap between them when plotted
        length = h.lengths(idx);
        C = squashed_rectangle(length/1.05,width*sum(h.lengths),0.2,points);
        
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

    %Generate half of 'an ellipse with height h and width cap_percent*w'
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