% creates a simple output plot using dynamic_error_chain

% setup
warning('off');
geometry.elements = {GroupElement.ROTARY, GroupElement.ROTARY,...
                     GroupElement.LINEAR, GroupElement.LINEAR,...
                     GroupElement.LINEAR, GroupElement.ROTARY};
geometry.types = {'j', 'e', 'j', 'e', 'l', 'j'};
syms theta x;
assume(x, 'real')
assume(theta, 'real')
% the below fixes length errors; are the error_chain defaults wrong?
configuration.errors = [theta x];
configuration.links = [1];
configuration.joints = [pi/2 0 0];
order = 2;
truncate_to = 2;
ranges = {-pi:pi/12:pi, 1/3};

% generate data
[joints, terms, sep_terms, errors, orders] = dynamic_error_chain(geometry,...
                                                                 configuration,...
                                                                 order,...
                                                                 truncate_to,...
                                                                 ranges);

% create baseline plot
clf;
hold on
line(joints(:,1), joints(:,2), 'Marker', 'o', 'Color', 'black');
xlim([-0.5 0.5]);
ylim([0 2]);
% plot error terms
legend_names = cell(1, length(errors)+1);
legend_names{1} = 'Nominal';
for t = 1:length(errors)
    plot(errors{t}(:,1), errors{t}(:,2));
    legend_names{t+1} = char(terms{t});
end
legend(legend_names);
hold off
display_terms(terms, sep_terms);

% fn to display existent terms
function display_terms(expected_terms, sep_error_terms)
    for term = 1:length(expected_terms)
        line = [char(expected_terms{term}), ': ', char(sep_error_terms{term}),';'];
        disp(line);
    end
end
