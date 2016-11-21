function update = check_sysplotter_dependencies(sys,shch)
% Determine which aspects of the system and path need to be recalculated

% Load the pathnames
configfile = './sysplotter_config';
pathnames = load(configfile);

% Verify that the system and the shape change both take two arguments
if ~isempty(sys) && nargin(sys) ~= 2;
    error('System function does not take two arguments')
end

if ~isempty(shch) && nargin(shch) ~= 2;
    error('Shape change function does not take two arguments')
end

	if ~isempty(sys)

		% Re-initialize system if necessary
		update.sys_init = depcheck([sys '.m'],pathnames);

		% Re-initialize shape change if necessary
		update.shch_init = ~isempty(shch) && depcheck([shch '.m'],pathnames);

		% Recalculate system details if necessary
		update.sys_calc = update.sys_init || depcheck('sys_calcsystem.m',sys);

		% Recalculate the path details if necessary
		update.shch_calc = (update.shch_init || update.sys_calc || depcheck('sys_calcpath.m',sys,shch));
	
	else
		
		[update.sys_init, update.shch_init,update.sys_calc,update.shch_calc] = deal(0);
		
	end
	
end