function varargout = absolute_feval(target, varargin)
%Evaluate the function whose mfile path is given by target, using varargin
%inputs

% Decide how many outputs we need
varargout = cell(nargout,1);

% Get the components of the target
[t_path, t_name, t_ext] = fileparts2(target);

% change to the target directory, run the function, then return to original
% directory
current_dir = pwd;
cd(t_path);

try % if this doesn't go through, will cause error below and drop back one level in the program
    filefunction = str2func(t_name);
	[varargout{:}] = filefunction(varargin{:});
catch ME
	cd(current_dir)
	rethrow(ME) % treat the error from feval as if it had happened, but first return to previous working directory
end


cd(current_dir);

end
