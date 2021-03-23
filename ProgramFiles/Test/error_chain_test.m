% test script for error_chain.m
% may also be used a reference for use

clear; clf;

% test error handling
disp('Errors:')
% no geometry
try
    error_chain();
catch e
    disp(e.message);
end
% empty geometry
geometry = struct;
try
    error_chain(geometry);
catch e
    disp(e.message);
end
% incorrect geometry elements
geometry.elements = {1 2 3 4};
geometry.types = cell(4);
try
    error_chain(geometry);
catch e
    disp(e.message);
end
% incorrect geometry types
geometry.elements = {GroupElement.LINEAR};
geometry.types = {1 2 3 4};
try
    error_chain(geometry);
catch e
    disp(e.message);
end

% construct error chains
disp([newline 'Rotary chain testing']);
geometry.elements = {GroupElement.ROTARY, GroupElement.ROTARY,...
                     GroupElement.LINEAR};
geometry.types = {'j', 'e', 'l'};
syms theta;
assume(theta, 'real');
configuration.errors = theta;
warning('on');
disp('Warnings about order, truncate_to:');
[~,~,expected_terms,sep_error_terms] = error_chain(geometry, configuration);
warning('off');
disp('First order error terms (default order, truncation):');
display_terms(expected_terms,sep_error_terms);
disp('First order error terms (correct order, truncation):');
[~,~,expected_terms,sep_error_terms] = error_chain(geometry, configuration, 2, 1);
display_terms(expected_terms,sep_error_terms);
disp('Second order error terms:');
[~,~,expected_terms,sep_error_terms] = error_chain(geometry, configuration, 3, 2);
display_terms(expected_terms,sep_error_terms);

disp([newline 'Rotary-Prismatic chain testing']);
geometry.elements = {GroupElement.ROTARY, GroupElement.ROTARY,...
                    GroupElement.LINEAR, GroupElement.LINEAR,...
                    GroupElement.LINEAR};
geometry.types = {'j','e','j','e','l'};
syms x;
assume(x, 'real');
configuration.errors = [theta x];
disp('First order error terms:');
[~,~,expected_terms,sep_error_terms] = error_chain(geometry, configuration, 2, 1);
display_terms(expected_terms,sep_error_terms);
disp('Second order error terms:');
[~,~,expected_terms,sep_error_terms] = error_chain(geometry, configuration, 3, 2);
display_terms(expected_terms,sep_error_terms);

% test N-link-chain outputs
disp([newline 'Default output testing']);
geometry.elements = {GroupElement.ROTARY, GroupElement.LINEAR,...
                     GroupElement.ROTARY, GroupElement.LINEAR,...
                     GroupElement.ROTARY};
geometry.types = {'j','l','j','l','j'};
configuration.joints = [pi/2, -pi/4, 0];
[~, ~, ~, ~, joint_configurations, J] = error_chain(geometry, configuration);
hold on
line(joint_configurations(:,1), joint_configurations(:,2), 'Color', 'black', 'Marker', 'o')
for j = 1:size(joint_configurations, 1)
    for c = 1:size(J, 2)
        quiver(joint_configurations(j,1), joint_configurations(j,2), J(1,c,j), J(2,c,j));
    end
end
xlim([0 3.5]);
ylim([0 3.5]);
hold off

% fn to display existent terms
function display_terms(expected_terms, sep_error_terms)
    for term = 1:length(expected_terms)
        line = [char(expected_terms{term}), ': ', char(sep_error_terms{term}),';'];
        disp(line);
    end
end

