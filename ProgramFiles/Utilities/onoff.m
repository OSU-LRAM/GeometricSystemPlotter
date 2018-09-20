function response = onoff(input)
% ONOFF returns logical 1 and 0 for strings that are 'on' and 'off' and
% vice versa (case-insensitive).  Also works for serial object statuses
% (i.e. 'open' or 'closed'), 'yes' and 'no', and dialog boxes ('ok' and
% 'cancel').

% error checking
narginchk(1, 1)

% if the input is a cell array of strings
if iscell(input)
    % recursively calls this function
    response = cellfun(@onoff, input, 'UniformOutput', false);

    % if the input was a cell array of strings, turn it into a matrix
    % (since it is now a cell array of logicals)
    if iscellstr(input)
        % converts it
        response = cell2mat(response);
    end

elseif ischar(input)
    % checks the input
    switch lower(input)
        case {'on', 'open', 'yes', 'ok'}
            response = true;

        case {'off', 'closed', 'no', 'cancel'}
            response = false;

        otherwise
            % errors
            error('Unknown string to convert.')
    end

elseif ~isempty(input) && ((isnumeric(input) && all(~isnan(input))) || islogical(input))
    % if it was a scalar
    if isscalar(input)
        if input
            response = 'on';
        else
            response = 'off';
        end

    else
        % recursively calls this function with the cell array form
        response = onoff(num2cell(input));
    end

else
    % error as nothing else can work
    error('onoff is defined for non-NaN non-empty arrays, strings and cell arrays only')
end