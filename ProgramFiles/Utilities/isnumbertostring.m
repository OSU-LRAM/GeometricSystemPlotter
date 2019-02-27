function [output] = isnumbertostring(value,n)
%isnumberself: This function will return [] if the input isnt a number and
%will return the string of that number if it is
% Inputs:
% Value: The number/string to be evaluated
% n: precision to be used if it is a number

if isnumeric(value)
    output = num2str(value,n);
else
    output = [];
end
end

