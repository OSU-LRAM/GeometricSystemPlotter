function update = check_sysplotter_dependencies(sys,shch)
% Determine which aspects of the system and path need to be recalculated

	if ~isempty(sys)

		% Re-initialize system if necessary
		update.sys_init = depcheck([sys '.m']);

		% Re-initialize shape change if necessary
		update.shch_init = ~isempty(shch) && depcheck([shch '.m']);

		% Recalculate system details if necessary
		update.sys_calc = update.sys_init || depcheck('sys_calcsystem.m',sys);

		% Recalculate the path details if necessary
		update.shch_calc = (update.shch_init || update.sys_calc || depcheck('sys_calcpath.m',sys,shch));
	
	else
		
		[update.sys_init, update.shch_init,update.sys_calc,update.shch_calc] = deal(0);
		
	end
	
end