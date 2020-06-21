function [frame_zero,J_zero,link_zero] = N_link_conversion_factors(C)
%%%%%%%
% This is a helper-function for N_link_chain. 
%
% N_link_chain takes in an argument that specifies which link on the chain
% (or another frame, such as a center-of-mass frame or one loaded from a
% system file) should be used as the base frame.
%
% All of the N_link_chain calculations are made with the lowest-numbered
% (leftmost, or 'tail') link as the base frame.
%
% Transforming the positions of the links relative to the base frame is
% straightforward: We simply multiply each transformation by the inverse of
% the transformation from the old base frame to the new base frame
%
% To transform the link Jacobians into the new coordinates, we need the
% Jacobian from joint velocities to body velocity of new frame, with
% first link fixed. This calculation depends on how the kinematic map
% from the original base frame to the new base frame is calculated.
%
% The code in this file handles the special cases required for finding the
% transformation and Jacobian from the old base frame to the new base frame
% for different specifications of the new base frame.
%
% The actual transformation of the link matrices and Jacobians is handled
% by N_link_conversion.

%%%%%%%%%%%%%%
% Get some basic information about the chain whose kinematics we're
% calculating

% Extract elements from structure array (some functions like size do not
% work properly if applied to struct contents without a breakout)
chain_m = C.chain_m;
jointchain_m = C.jointchain_m;
links_m = C.links_m;
%joints_m = C.joints_m;
links_v = C.links_v;
joints_v = C.joints_v;
%jointangles = C.jointangles;
jointangles_c = C.jointangles_c;
linklengths = C.linklengths;
shapeparams = C.shapeparams;
modes = C.modes;
J_temp = C.J_temp;
baseframe = C.baseframe;


% Number of links and joints
N_links = size(chain_m,3);
M_joints = N_links-1;

% Decide if there is an even or odd number of joints
N_odd = mod(N_links,2);

% Total length for use in taking weighted averages
L = sum(linklengths);

%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%
% If baseframe is specified as a non-cell object, wrap it into a cell array
if ~iscell(baseframe)
    baseframe = {baseframe};
end


% Depending on which baseframe option is specified, use different methods
% to calculate transform and Jacobian to use in the conversion
for idx_baseframe = 1:numel(baseframe)
    
    if ischar(baseframe{idx_baseframe}) || isscalar(baseframe{idx_baseframe})

        switch baseframe{idx_baseframe}

            % Places the reference frame at the midpoint of the lowest-index link on the chain
            case 'tail'

                % Identify the first link
                link_zero = 1;

                % Extract the transform from the end link to itself 
                frame_zero = eye(3);       

                % Jacobian for base link is zero
                J_zero = zeros(size(J_temp{1}));

            % Places the reference frame at the *end* of the lowest-index link on the chain
            case {'tail-tip','tail-end'}

                % Identify the end link
                link_zero = 1;

                % The transform from the midpoint of the first link to its end is
                % the inverse of the half-link transform
                frame_zero = inv(links_m(:,:,link_zero));

                %%%%%%%%
                % Jacobian to new frame is Jacobian of first link, but with an
                % adjoint-inverse transform by the half-link to get to the end

                halfstep = Adjinv(frame_zero);
                J_zero = halfstep * J_temp{1};


            % Places the reference frame at the middle of the chain, splitting the
            % difference between the two middle links if there is an even number of
            % links
            case {'centered','center','midpoint-tangent'}


                % If there is an odd number of links, the middle link is the center of the
                % chain
                if N_odd

                    % Identify the middle link
                    link_zero = ceil(N_links/2);

                    % Extract the transform from the end link to the middle link
                    frame_zero = chain_m(:,:,link_zero);

                    %%%%%%
                    % Jacobian for conversion is Jacobian of link zero
                    J_zero = J_temp{link_zero};

                % If there is an even number of links, the midpoint is at the middle joint,
                % rotated by half of that joint's angle
                else

                    % Identify the middle joint
                    joint_zero = ceil(M_joints/2);

                    % Extract the position of the middle joint, and multiply it by half of
                    % the transform associated with its joint angle
                    frame_zero = jointchain_m(:,:,joint_zero) * ...
                                        vec_to_mat_SE2(joints_v(joint_zero,:)/2);

                    %%%%%%%%
                    % Jacobian for conversion is Jacobian of link before it, with a
                    % half-step to get to the end of the link, and a half rotation
                    % to average orientation between the two links
                    halfstep = Adjinv(links_v(joint_zero,:));
                    halfrotation = Adjinv(joints_v(joint_zero,:)/2);
                    J_zero = halfrotation * halfstep * J_temp{joint_zero};

                    % Account for centered frame being half-sensitive to middle joint
                    J_zero = J_zero + .5*[zeros(2,size(modes,2));modes(joint_zero,:)];

                    %J_zero(:,joint_zero) = [0; 0; .5]; 


                end





            % Places the reference frame at the midpoint of the highest-numbered link in the chain
            case 'head'

                % Identify the end link
                link_zero = N_links;

                % Extract the transform from the end link to the end link
                frame_zero = chain_m(:,:,link_zero);

                %%%%%%%%
                % Jacobian to new frame is Jacobian of last link
                J_zero = J_temp{end};

            % Places the reference frame at the *end of* the highest-numbered link
            % in the chain
            case 'head-tip'

                % Identify the end link
                link_zero = N_links;

                % Extract the transform from the first link to the end link, and
                % multiply it by the half-link transformation on that link
                frame_zero = chain_m(:,:,link_zero)*links_m(:,:,link_zero);

                %%%%%%%%
                % Jacobian to new frame is Jacobian of last link, but with an
                % adjoint-inverse transform by the half-link to get to the end
                halfstep = Adjinv(links_v(end,:));
                J_zero = halfstep * J_temp{end};

            % Places the reference frame at center of mass and average orientation
            % of the links, using link-lengths as weighting terms
            case 'com-mean'

                % Convert link positions to row form
                chain = mat_to_vec_SE2(chain_m);

                % Substitute a cumulative sum of the joint angles for the
                % orientations of the links: We're not looking for the mean of the
                % SO(2)-modulated orientations, we're looking for the mean angular
                % displacement from the base link.
                chain(:,3) = jointangles_c;

                % Take a weighted average of the link positions
                CoM = sum(diag(linklengths)*chain)/L;


                % Place the new frame at this location
                frame_zero = vec_to_mat_SE2(CoM);

                %%%%%%%%%%%
                % The Jacobian of the weighted average of frames is the
                % weighted average of their Jacobian (by the commutativity of
                % sumation and derivation operations).

                % Multiply each link's Jacobian by its link length
                J_weighted = J_temp;
                for idx = 1:numel(J_weighted)
                    J_weighted{idx} = TeLg(chain(idx,:)) * J_temp{idx} * linklengths(idx);
                end

                % Sum the weighted Jacobians
                J_zero = J_weighted{1};
                for idx = 2:numel(J_weighted)
                    J_zero = J_zero + J_weighted{idx};
                end
        %        J_zero = sum(cat(3,J_weighted{:}),3);

                % Divide by the total length to g
                J_zero = J_zero/L;  

                % Bring into local coordinates
                J_zero = TgLginv(frame_zero)*J_zero;


            %%%%%%%%%%%
            % These are cases that can be used as modifiers on other baseframes

            % Places the reference frame at the proximal end of the link to
            % which it is attached
            case 'start'

                % Make sure that a link was previously specified
                if ~exist('link_zero','var')
                    error('Baseframe string ''start'' may only be used as a secondary baseframe specification, where the primary specification attaches the baseframe to a link.');
                end

                % The transform from the midpoint of the first link to its proximal end is
                % the inverse of the half-link transform
                frame_step = inv(links_m(:,:,link_zero));

                % compose the frame step with the previously-computed
                % frame_zero
                frame_zero = frame_zero*frame_step; %#ok<MINV>

                %%%%%%%%
                % Jacobian to new frame is previously-calculated Jacobian, but with an
                % adjoint-inverse transform by the half-link to get to the end

                halfstep = Adjinv(frame_step);
                J_zero = halfstep * J_zero;

            % Places the reference frame at the distal end of the link to which
            % it is attached
            case 'end'

                % Make sure that a link was previously specified
                if ~exist('link_zero','var')
                    error('Baseframe string ''end'' may only be used as a secondary baseframe specification, where the primary specification attaches the baseframe to a link.');
                end

                % Transform from the midpoint to the distal end of the link
                frame_step = links_m(:,:,link_zero);

                % compose the frame step with the previously-computed
                % frame_zero
                frame_zero = frame_zero*frame_step;


                %%%%%%%%
                % Jacobian to new frame is previously-calculated Jacobian, but with an
                % adjoint-inverse transform by the half-link to get to the end

                halfstep = Adjinv(frame_step);
                J_zero = halfstep * J_zero;
            
            
            
        %%%%%%%%%%%
        % Several possible cases are handled within this "otherwise" case,
        % because they need if/else logic
        otherwise

            %%%%%%%
            % If a number is provided as the input, place the baseframe on that
            % link
            % Check for numeric baseframe specification
            if isnumeric(baseframe{idx_baseframe}) && isscalar(baseframe{idx_baseframe})

                % Make sure that base frame specification is actually in
                % the range of valid link numbers
                if baseframe{idx_baseframe} <= N_links &&  baseframe{idx_baseframe} > 0

                    % Set the base frame as numerically specified in the input
                    link_zero = baseframe{idx_baseframe};

                    % Extract the transform from the end link to the end link, and
                    % multiply it by the half-link transformation on that link
                    frame_zero = chain_m(:,:,link_zero);

                    %%%%%%%%
                    % Jacobian to new frame is Jacobian of nth link
                    J_zero = J_temp{baseframe{idx_baseframe}};

                else

                    error (['Baseframe specification ' num2str(baseframe{idx_baseframe}) ' is not in the range of link numbers'])

                end
                

            
            
            %%%%%%%%
            % If the provided string is the name of a sysf_ system file, pull
            % the transformation to minimum-perturbation coordinates from that file            
            else

                %%%%%
                % Robustly check if the baseframe string is a sysf_ file that
                % has had minimum-perturbation coordinates calculated

                %%%
                % First, strip off sysf_ if it was included by the user

                if strncmp(baseframe{idx_baseframe},'sysf_',5)
                    baseframe1 = baseframe{idx_baseframe}(6:end);
                else
                    baseframe1 = baseframe{idx_baseframe};
                end

                %%%
                % Second, put sysf_ onto the front of the string

                baseframe1 = ['sysf_' baseframe1]; %#ok<AGROW>

                %%%
                % Third, check if there is a sysf_ file with this name

                load('sysplotter_config.mat','inputpath')      % Get the path to the current UserFiles folder

                sysf_filepath = fullfile(inputpath, 'Systems',...  % Build a string to a sysf_ file with the baseframe name
                    [baseframe1 '.m']);

                sysf_file_present =  exist(sysf_filepath, 'file'); % Test for this file being present

                %%%
                % Fourth, if there is a sysf_ file of the appropriate name, see
                % if there is an associated file with the results of the system
                % calculations
                if sysf_file_present

                    % Pulling the minimum-perturbation coordinates does not
                    % work with symbolic input, because it requires an
                    % interpolation of an array at the input points. (Given
                    % that minimum-perturbation coordinates are solved for
                    % numerically, not much is lost by not being able to
                    % symbolically rebase to that frame
                    if isa(chain_m,'sym')
                        error('Rebasing the chain to minimum-perturbation coordinates is not supported for symbolic input');
                    end

                    % Build a string to the filename where calculations from
                    % this system are stored
                    calcfilePath = fullfile(inputpath, 'sysplotter_data',...
                        [baseframe1 '_calc.mat']);

                   % Check that the calculated system file exists
                   if exist(calcfilePath, 'file') 

                       % Import the datafile of the specified system
                        load(calcfilePath,'s')


                        %%%%
                        % Get the Jacobian from the baseframe used by the
                        % system to the minimum-perturbation baseframe



                        %%%
                        % Interpolate the specified shape into the grids to
                        % extract the transform from the system's baseframe to
                        % its minimum-perturbation baseframe and the Jacobian
                        % of this transform

                        frame_mp = zeros(3,1);              % x y theta values of transformation
                        J_mp = zeros(3,numel(shapeparams)); % x y theta rows and joint angle columns of Jacobian

                        shapeparams_cell = num2cell(shapeparams); % put joint angles into a cell array

                        % Iterate over x y theta components
                        for idx = 1:numel(frame_mp)

                            % Interpolate the shape angles into the grid and
                            % extract the corresponding component of the
                            % transformation to m-p coordinates
                            frame_mp(idx) = interpn(s.grid.eval{:},...                  % system evaluation grid
                                                    s.B_optimized.eval.Beta{idx},...    % Components of the transfomation to m-p coordinates
                                                    shapeparams_cell{:});               % Current shape

                            % Scale the x and y components of the frame
                            % location 
                            if any(idx==[1 2])
                                frame_mp(idx) = L*frame_mp(idx)/s.geometry.length;
                            end

                            % Iterate over shape components for Jacobian
                            for idx2 = 1:numel(shapeparams)


                                % Interpolate the shape angles into the grid
                                % and extract the corresponding component of
                                % the transformation to m-p coordinates
                                J_mp(idx,idx2) = interpn(s.grid.eval{:},...                  % system evaluation grid
                                                        s.B_optimized.eval.gradBeta{idx,idx2},...    % Components of the transfomation to m-p coordinates
                                                        shapeparams_cell{:});               % Current shape

                                % Scale the x and y components of the frame
                                % Jacobian
                                if any(idx==[1 2])
                                    J_mp(idx,idx2) = L*J_mp(idx,idx2)/s.geometry.length;
                                end
                            end



                        end



                        % Convert the transforms that were just found into an
                        % SE(2) matrix and a Jacobian for the body velocity
                        % of the new frame relative to the frame in which it
                        % is defined
                        frame_mp = vec_to_mat_SE2(frame_mp);
                        J_mp = TgLginv(frame_mp) * J_mp;
                        
                        
                        % Concatenate this transform and Jacobian with any
                        % previously-calculated values
                        frame_zero = frame_zero * frame_mp;
                        J_zero = (Adjinv(frame_mp) * J_zero) + J_mp;

%                        %%%%
%                         % Get the conversion factors from the tail link to the
%                         % baseframe used by the system
% 
%                         % Get the geometry specification in this file
%                         if isfield(s,'geometry')
%                             if isfield(s.geometry,'baseframe')
%                                 baseframe_original = s.geometry.baseframe;
%                             else
%                                 error('Baseframe is not specified in system geometry structure')
%                             end
%                         else
%                             error('System does not have a geometry structure')
%                         end
% 

%                         % Get conversion factors for the original baseframe
%                         [frame_original,J_original] =...
%                             N_link_conversion_factors(C);
% %                         chain_m,...
% %                                                      jointchain_m,...
% %                                                      links_m,...
% %                                                      joints_m,...
% %                                                      links_v,...
% %                                                      joints_v,...
% %                                                      jointangles,...
% %                                                      linklengths,...
% %                                                      shapeparams,...
% %                                                      modes,...
% %                                                      J_temp,...
% %                                                      baseframe_original);
% 
%                         % Combine original conversion with minimum-perturbation
%                         % conversion
% 
%                         % Frame construction is straightforward concatenation
%                         % of SE(2) transforms
%                         frame_zero = frame_original * frame_mp; 
% 
%                         % Combination of body-velocity Jacobians is by taking
%                         % the adjoint-inverse of the proximal Jacobian by the
%                         % distal transformation, then adding the distal
%                         % Jacobian
%                         J_zero = (Adjinv(frame_mp) * J_original) + J_mp;

                   end

                else
                    error (['Baseframe specification ' baseframe{idx_baseframe} ' is not a valid option, and there is no system file with that name'])
                end
                
            end
        end
            
        %%%%%%%
        % If a transformation matrix is provided as the input, apply it
        % to the base frame (stacked onto any previous base frame)

    elseif isnumeric(baseframe{idx_baseframe})


        if ~exist('J_zero','var')

            % Identify the first link
            link_zero = 1;

            % Extract the transform from the end link to itself 
            frame_zero = eye(3);       

            % Jacobian for base link is zero
            J_zero = zeros(size(J_temp{1}));

        end

        % Step by the provided transformation
        frame_step = baseframe{idx_baseframe};

        % Apply the step relatve to the pre-existing frame zero
        frame_zero = frame_zero * frame_step;

        % Use the adjoint-inverse of the step to modify the
        % Jacobian
        J_zero = Adjinv(frame_step) * J_zero;

    end
            
                     
end
      
    % Update J_temp to be the J_zero computed at this stage of the loop
%    J_temp = J_zero;
    
if ~exist('link_zero','var')
    link_zero = 0;
end

end


    