function sys_update(sys,shch,stretch,progress,handles)
    %Make sure data in the plot file is up to date
	
	% declare the data directory
	configfile = './sysplotter_config';
	load(configfile);
    
    
	if exist('progress','var')
	
		%Progress bar stages, step size, and state
		prog_stages = 5;
		prog_step = 1/prog_stages;
		prog_state = 0;

		%Advance the progress bar
		prog_state = prog_state + 1;
		waitbar2a(prog_state*prog_step,progress,'Checking Dependencies')
	
	end
	
	%convert empty matrix input to 'null' string
	if isempty(shch)

        shch = 'null';

	end

	% Determine which components need to be updated
	if exist('handles','var')
		[update] = refresh_runinfo_Callback([], [], handles);
	else
		[update] = decide_components(sys,shch);
	end
	
    %%%%%%%%%%
	% Recalculate the details if necessary

    % Load a structure with the path names
    pathnames = load(configfile);
	
	if update.sys_init

        s = absolute_feval(fullfile(syspath, sys),'initialize',pathnames); %#ok<NASGU>
        
		save(fullfile(datapath, sys),'s')
	end
	
	if update.shch_init

       p = absolute_feval(fullfile(shchpath, shch),'initialize',pathnames); %#ok<NASGU>
        
		save(fullfile(datapath, shch),'p')
	end
		
	if exist('progress','var')
		%Advance the progress bar
		prog_state = prog_state + 1;
		waitbar2a(prog_state*prog_step,progress,'Updating System')
	end
	
	if update.sys_calc
		sys_calcsystem('calculate',sys);
    end
	
	if exist('progress','var')
		%Advance the progress bar
		prog_state = prog_state + 1;
		waitbar2a(prog_state*prog_step,progress,'Updating Path')
	end
	
	if update.shch_calc
		sys_calcpath('calculate',sys,shch);
	end
	



	if exist('progress','var')
		%Advance the progress bar
		prog_state = prog_state + 1;
		waitbar2a(prog_state*prog_step,progress,'Plotting Data')
	end

end

    




