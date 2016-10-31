function rel = isrel(str)
% Decide if the input str is a relative or absolute path, by comparing
% first four characters with those from pwd

teststr = pwd;

rel = ~strncmp(str,teststr,4);

end