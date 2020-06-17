function [frame_zero,J_zero] = backbone_conversion_factors(h,...
                                                           J,...
                                                           shapeparams,...
                                                           baseframe,...
                                                           L)
%%%%%%%
% This is a helper-function for backbone. 
%
% backbone takes in an argument that specifies which point on the backbone
% (or another frame, such as a center-of-mass frame or one loaded from a
% system file) should be used as the base frame.
%
% All of the backbone calculations are made with the midpoint
% (s=0) frame as the base frame.
%
% Transforming the positions of the backbone pieces relative to the base
% frame is straightforward: We simply multiply the transformation function
% for the backbone pieces by the inverse of the transformation from the old
% base frame to the new base frame
%
% To transform the link Jacobians into the new coordinates, we need the
% Jacobian from shape velocities to body velocity of new frame, with
% the center frame fixed. This calculation depends on how the kinematic map
% from the original base frame to the new base frame is calculated.
%
% The code in this file handles the special cases required for finding the
% transformation and Jacobian from the old base frame to the new base frame
% for different specifications of the new base frame.
%
% The actual transformation of the link matrices and Jacobians is handled
% by backbone_conversion.

% Only calculate Jacobians if requested in the output
if nargout > 1
    calc_J = 1;
else
    calc_J = 0;
end
    
%%%%%%%%%%%%%%%%%%%
% If baseframe is specified as a non-cell object, wrap it into a cell array
if ~iscell(baseframe)
    baseframe = {baseframe};
end

% Depending on which baseframe option is specified, use different methods
% to calculate transform and Jacobian to use in the conversion
% Depending on which baseframe option is specified, use different methods
% to calculate transform and Jacobian to use in the conversion
for idx_baseframe = 1:numel(baseframe)
    
    if ischar(baseframe{idx_baseframe}) || isscalar(baseframe{idx_baseframe})
        
        switch baseframe{idx_baseframe}
    
            % Places the reference frame at the leftmost (lowest s) point on the
            % backbone
            case 'tail'

                % Identify the leftmost portion of the body
                s_zero = -0.5;

                % Extract the transform for this frame 
                frame_zero = vec_to_mat_SE2(h(s_zero));       

                % Jacobian for tail frame
                if calc_J
                    J_zero = TgLginv(frame_zero)*J(s_zero);
                end

            % Places the reference frame at the middle of the chain, splitting the
            % difference between the two middle links if there is an even number of
            % links
            case {'centered','center','midpoint-tangent'}


                % Identify the middle portion of the body
                s_zero = 0;

                % Extract the transform for this frame 
                frame_zero = vec_to_mat_SE2(h(s_zero));       

                % Jacobian for middle frame
                if calc_J
                    J_zero = TgLginv(frame_zero)*J(s_zero);
                end

            % Places the reference frame on the rightmost (highest s) point on the
            % backbone
            case {'head','head-tip'}

                % Identify the leftmost portion of the body
                s_zero = 0.5;

                % Extract the transform for this frame 
                frame_zero = vec_to_mat_SE2(h(s_zero));       

                % Jacobian for tail frame
                if calc_J
                    J_zero = TgLginv(frame_zero)*J(s_zero);
                end


            % Places the reference frame at center of mass and average orientation
            % of the backbone, assuming constant mass-per-s density
            case 'com-mean'

                % Sample the backbone at a dense set of points
                h_points = h(linspace(-.5,.5,100));

                % Average the x, y, and theta components. The theta value here is
                % integrated local rotation, so we don't need to worry about
                % unwrapping it
                CoM = mean(h_points,2);

                % Place the new frame at this location
                frame_zero = vec_to_mat_SE2(CoM);

                %%%%%%%%%%%
                % The Jacobian of the weighted average of frames is the
                % weighted average of their Jacobian (by the commutativity of
                % sumation and derivation operations).

                if calc_J

                    % Sample the Jacobian at a dense set of points
                    J_points = J(linspace(-.5,.5,100));

                    % Average the Jacobians
                    J_zero = mean(J_points,3);

                    % Bring into local coordinates
                    J_zero = TgLginv(frame_zero)*J_zero;

                end        

            %%%%%%%%%%%
            % Several possible cases are handled within this "otherwise" case,
            % because they need if/else logic
            otherwise

                %%%%%%%
                % If a number is provided as the input, place the baseframe on that
                % link
                % Check for numeric baseframe specification
                if isnumeric(baseframe) && isscalar(baseframe{idx_baseframe})

                    % Make sure that base frame specification is actually in
                    % the range of valid link numbers
                    if baseframe >= -0.5 &&  baseframe <= 0.5

                        % Set the base frame as numerically specified in the input
                        s_zero = baseframe;

                        % Extract the transform for this frame 
                        frame_zero = vec_to_mat_SE2(h(s_zero));       

                        % Jacobian for middle frame
                        if calc_J
                            J_zero = TgLginv(frame_zero)*J(s_zero);
                        end

                    else

                        error (['Baseframe specification ' num2str(baseframe) ' is not in the range of -0.5 to 0.5'])

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

                            if calc_J
                                J_mp = zeros(3,numel(shapeparams)); % x y theta rows and joint angle columns of Jacobian
                            end

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

                                if calc_J
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

                            end

                            % Convert the transforms that were just found into an
                            % SE(2) matrix and a body-velocity Jacobian
                            frame_mp = vec_to_mat_SE2(frame_mp);

                            if calc_J
                                J_mp = TgLginv(frame_mp) * J_mp;
                            end

                            %%%%
                            % Get the conversion factors from the tail link to the
                            % baseframe used by the system

                            % Get the geometry specification in this file
                            if isfield(s,'geometry')
                                if isfield(s.geometry,'baseframe')
                                    baseframe_original = s.geometry.baseframe;
                                else
                                    error('Baseframe is not specified in system geometry structure')
                                end
                            else
                                error('System does not have a geometry structure')
                            end


                            % Get conversion factors for the original baseframe
                            if calc_J
                                [frame_original,J_original] =...
                                    backbone_conversion_factors(h,J,shapeparams,baseframe_original,L);
                            else
                                frame_original =...
                                    backbone_conversion_factors(h,J,shapeparams,baseframe_original,L);                        
                            end

                            % Combine original conversion with minimum-perturbation
                            % conversion

                            % Frame construction is straightforward concatenation
                            % of SE(2) transforms
                            frame_zero = frame_original * frame_mp; 

                            % Combination of body-velocity Jacobians is by taking
                            % the adjoint-inverse of the proximal Jacobian by the
                            % distal transformation, then adding the distal
                            % Jacobian
                            if calc_J
                                J_zero = (Adjinv(frame_mp) * J_original) + J_mp;
                            end

                       end

                    else
                        error (['Baseframe specification ' baseframe ' is not a valid option, and there is no system file with that name'])
                    end

                end

        end
        
        elseif isnumeric(baseframe{idx_baseframe})


        if ~exist('J_zero','var')

            % Extract the transform from the end to itself 
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

end
    