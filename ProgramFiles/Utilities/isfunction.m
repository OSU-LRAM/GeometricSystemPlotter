function truefalse = isfunction(inputstr)

try
    nargin(inputstr);
    truefalse = true;
catch
    truefalse = false;
end