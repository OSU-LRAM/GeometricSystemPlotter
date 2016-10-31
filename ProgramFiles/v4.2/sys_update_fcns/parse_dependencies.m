function [depstructure, full_target] = parse_dependencies(target,varargin)
% For files that may have dependencies specified by local paths, or default
% outputs, make them into absolute paths

% load the list of paths
configfile = './sysplotter_config';
load(configfile)
	

% Generate base paths for local dependencies based on the file type, along
% with default products to insert if no product is given
if strncmp('sysf_',target,5)
	
	dep_base = syspath;
	default_product = {fullfile(datapath,[target(1:end-2) '.mat'])};
		
elseif strncmp('shchf_',target,6);

	dep_base = shchpath;
	default_product = {fullfile(datapath,[target(1:end-2) '.mat'])};

elseif strncmp('stretchf_',target,9);
	
	dep_base = stretchpath;
	default_product = {};

else
	
	dep_base = pwd;
	default_product = {};

end

% Form the full target
full_target = make_path_absolute(target,dep_base);

% Get the components of the target
[t_path, t_name, t_ext] = fileparts2(full_target);

% Get the raw dependency structure for the system
current_dir = pwd;
cd(t_path);
try % if this doesn't go through, will cause error below and drop back one level in the program
	depstructure = feval(t_name,'dependency',varargin{:});
catch
%	warning('Evaluation of absolutely defined dependency failed')
end
cd(current_dir);

% Turn all paths in dependencies into absolute paths
depstructure.dependency = cellfun(@(str) make_path_absolute(str,dep_base),depstructure.dependency,'UniformOutput',false);

% Give the structure a product if it doesn't have one (note that an empty
% product is the null product, while a missing field should be filled in)
if ~isfield(depstructure,'product')
	depstructure.product = default_product;
end

end