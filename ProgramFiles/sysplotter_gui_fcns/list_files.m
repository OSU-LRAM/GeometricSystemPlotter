function [displaynames, filenames] = list_files(dirname,prefix)
% List the systems or path files valid for use with the gui
%
% Loop over the contents of the target directory, asking each for its
% displayname


	% Make sure the Systems directory is in the path
	addpath(dirname)
	
	% Get a list of all system m-files
	D = dir([dirname '/' prefix '_*.m']);
	
	% Extract the names of the files
	[displaynames, filenames] = deal(cell(length(D),1));
	
	for i = 1:length(D)
		
        filefunction = str2func(D(i).name(1:end-2)); % filename into a function
        
		displaynames{i} = filefunction('name'); % Get system files to declare their names
		
		filenames{i} = D(i).name(1:end-2); % Record the filename of the system file
		
	end
	
	% Sort the names alphabetically by the display names
	[displaynames, nameI] = sort(displaynames);
	filenames = filenames(nameI);
	

end