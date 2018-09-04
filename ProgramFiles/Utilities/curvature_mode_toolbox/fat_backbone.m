function [B,h,J,J_full,frame_zero,J_zero] = fat_backbone(geometry,shapeparams,display)
% Generate a rounded-end thick curve corresponding to a continuous backbone
% shape
% 
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields used by this function are:
%
%       type: The kind of backbone function provided here. Used by backbone
%           to calculate the geometry of the backbone
%
%       function: A handle to a curvd_ function describing the system's
%           curvature function, or a cell vector of function handles to the
%           curvature modes
%
%       baseframe (optional): Specification of which point on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the backbone. Options for this field are
%           specified in the documentation of backbone.
%
%
%       length (optional): Total length of the backbone. If specified, the
%           elements of will be scaled such that their sum is equal to L.
%           If this field is not provided or is entered as an empty matrix,
%           then the backbone is taken as being of unit length.
%
%   shapeparams: A vector of the  input parameters taken by the curvature
%       function
%
%   display: Structure containing information about how the system should
%       be illustrated. Field required by this function is:
%
%       aspect: ratio of width to length of the system
%
%       sharpness: fraction of the backbone taken up by rounded tip
%
%
%
% Outputs:
%
%   B : Locus of points for the fat backbone. First two rows are x and y
%       values, third row is all ones (so that the columns can be
%       transformed by SE(2) matrices.
%
%   h : Function handle for h(s) returning location of point s on the
%       backbone see the function backbone for more information
%
%   J : Jacobians from joint angle velocities to velocity of each link
%           relative to  the base frame (the standard derivative of h).
%           These Jacobians are stored in a cell array, with one matrix per
%           cell.
%
%   J_full : Full Jacobians from body velocity of the base frame and joint
%           angle velocities to body velocity of each link relative to a
%           non-moving frame.
%
%   frame_zero : The transformation from the first link to the selected
%           baseframe
%
%   J_zero : The Jacobian from shape velocity to body velocity of the
%       selected baseframe, with the first link held fixed.


%%%%%%%%%%%
    % Only calculate Jacobians if requested in the output
    if nargout > 1
        calc_J = 1;
    else
        calc_J = 0;
    end


    % Fill in default display parameters if not provided
    if ~exist('display','var') || ~isfield(display,'aspect_ratio')
        display.aspect_ratio = 0.03;
    end
    if ~isfield(display,'sharpness')
        display.sharpness = 0.03; % This matches some shapes originally drawn in OmniGraffle
    end
    
    % If no length is specified, use unit length
    if ~isfield(geometry,'length') || isempty(geometry.length)
        geometry.length = 1;
    end
    
    % Request the function for the backbone, and whichever other functions
    % have been requested as passthroughs
    if calc_J
        [h,J,J_full] = backbone(geometry,shapeparams);
    else
        h = backbone(geometry,shapeparams);
    end
    
    % Turn this backbone function into a locus of points
    B = fatbackbone_locus(geometry.length,...
                          display.aspect_ratio,...
                          display.sharpness,...
                          h);
    
    

end






function footprint = fatbackbone_locus(L,aspect_ratio,sharpness,backbone_fun)

    % Figure out how much of the system is taken up by the end caps
	a = L*sharpness;  
	b = L*aspect_ratio/2;

    % Build the end caps as partial ellipses
	poscaptheta = linspace(pi/2,-pi/2,30)';
	poscapprim = [a*cos(poscaptheta) b*sin(poscaptheta)];
	
	negcaptheta = linspace(3*pi/2,pi/2,30)';
	negcapprim = [a*cos(negcaptheta) b*sin(negcaptheta)];

    % get the positions and orientations of points along the core
	corepoints = linspace(-0.5,0.5,200);
    core = backbone_fun(corepoints);
		
	% find the location of points along the perimeter
	posedge = zeros(length(corepoints),2);
	negedge = posedge;
	

	
	for i = 1:length(corepoints)
        
        transform = vec_to_mat_SE2(core(:,i));
		
		
		edgepoints = transform*[0 0;b -b;1 1];
		
		posedge(i,:) = edgepoints(1:2,1)';
		negedge(end+1-i,:) = edgepoints(1:2,2)';
		
		if i == 1
			negcap = (transform*[negcapprim';ones(1,length(negcaptheta))])';
			negcap(:,3) = [];
		end
		
		if i == length(corepoints)
			poscap = (transform*[poscapprim';ones(1,length(poscaptheta))])';
			poscap(:,3) = [];
		end

		
	end	
	
	

	footprint = [posedge;poscap;negedge;negcap]';
	
    % Augment this footprint with a row of ones so that it can be operated
    % on via SE(2) matrices
    footprint = [footprint;ones(1,size(footprint,2))];

end