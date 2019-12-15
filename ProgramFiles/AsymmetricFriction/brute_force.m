% TODO: find out how to un-hardcode this
%home = "/home/qkonyn/Programming_Projects/GeometricSystemPlotter/";
%addpath(home + "UserFiles/GenericUser/Systems");
%addpath(home + "ProgramFiles/Geometry");

s = sysf_two_link_lowRe;

friction_direction = 1; % todo: vary this
shape_step = 0.05;
% shape is alpha; shapechange is alpha dot
a = s.grid_range(1):shape_step:s.grid_range(2);
adot = -0.5:0.01:0.5;

% storage of results for each link
agree = cell(size(s.geometry.linklengths));
% how do I do this without hardcoding?
agree{1} = zeros(length(a), length(adot));
agree{2} = zeros(length(a), length(adot));

for i = 1:length(a)
    shape = a(i);
    
    % these two don't depend on alpha dot:
    % get A(alpha)
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape);
    % get Jfull(alpha)
    [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
    
    for j = 1:length(adot)
        shapechange = adot(j);
        
        % gcirc right
        body_vel = A * shapechange;

        for link = 1:length(s.geometry.linklengths)
            body_vel = J_full{link} * [body_vel; shapechange];
            dx_direction = sign(body_vel(1)); %index 1 is x, right?
            if(dx_direction == 0) % non-if way to do this?
                agreement = 1;
            else
                agreement = friction_direction * dx_direction;
            end
            agree{link}(i,j) = agreement;
        end
    end
end

%surface(a, adot, agree{1}); % link 1
%surface(a, adot, agree{2}); % link 2

figure('Name','Consistency Summary: Forward');
% summary: link 1 is the 1's digit in binary. link 2 is the 2's digit.
% so 00 would be agreeing with neither, 10 (= 2) would be agreeing with
% just link 2, 01 (= 1) would be agreeing with just link 1, and 11 (= 3)
% would be agreeing with both
summary = (agree{1} + 1)/2 + (agree{2} + 1);
surface(a, adot, summary);

% todo: DRY
friction_direction = -1;
s.physics.drag_coefficient = s.physics.drag_coefficient2;
shape_step = 0.05;
% shape is alpha; shapechange is alpha dot
a = s.grid_range(1):shape_step:s.grid_range(2);
adot = -0.5:0.01:0.5;

% storage of results for each link
agree = cell(size(s.geometry.linklengths));
% how do I do this without hardcoding?
agree{1} = zeros(length(a), length(adot));
agree{2} = zeros(length(a), length(adot));

for i = 1:length(a)
    shape = a(i);
    
    % these two don't depend on alpha dot:
    % get A(alpha)
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape);
    % get Jfull(alpha)
    [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
    
    for j = 1:length(adot)
        shapechange = adot(j);
        
        % gcirc right
        body_vel = A * shapechange;

        for link = 1:length(s.geometry.linklengths)
            body_vel = J_full{link} * [body_vel; shapechange];
            dx_direction = sign(body_vel(1)); %index 1 is x, right?
            if(dx_direction == 0) % non-if way to do this?
                agreement = 1;
            else
                agreement = friction_direction * dx_direction;
            end
            agree{link}(i,j) = agreement;
        end
    end
end

figure('Name','Consistency Summary: Backward');
summary = (agree{1} + 1)/2 + (agree{2} + 1);
surface(a, adot, summary);
