function [update, dep_update] = decide_components(sys,shch,handles)
% Decide which components will be rerun, based on the dependencies and the
% items selected in the gui

	% First, check the file dependencies
	dep_update = check_sysplotter_dependencies(sys,shch);
	
	
	if exist('handles','var')
		% Second, check the gui selection
		gui_update = collect_recalc_button_values(handles);

		% Modify the dependency updates with information from the gui
		update = dep_update;

		components = fieldnames(update);


		switch gui_update(end)

			% If the override is set to "auto", use the unmodified dependency check
			case 1

				% do nothing

			% If the override is set to "force", update everything
			case 2

				for i = 1:length(components)

					update.(components{i}) = 1;

				end

			% If the override is set to "block", don't update anything
			case 3

				for i = 1:length(components)

					update.(components{i}) = 0;

				end

			% If the override is set to "Individual", use the setting for each
			% component
			case 4

				relevance_list = {1 2 [1 3] [2 3 4]}; % which gui blocks are relevant to each execution stage
				
				for i = 1:length(components)

					update.(components{i}) = ...
						(dep_update.(components{i}) || any(gui_update(relevance_list{i}) == 2)) ...
						&& (gui_update(i) ~= 3); %dependency or gui calls for calculation, and does not block it

				end
		end
		
	else % If this function was not called from the gui, use only the dependencies
		
		update = dep_update;
		
	end
	
	% if the shape change is null, this beats even "force"
	if strcmp('null',shch)
		
		update.shch_init = 0;
		
	end

end