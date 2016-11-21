function str = make_path_absolute(str,base)
% If the path in str is a relative path, make it an absolute path starting
% from base

if isrel(str)
	str = fullfile(base, str);
end

end