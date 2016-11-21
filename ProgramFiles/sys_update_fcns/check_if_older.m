function older = check_if_older(target,reference)
% Check if the target is older than the reference

	% Get the info about the target and reference
	target_info = batch_dir(target);
	reference_info = batch_dir(reference);

	% Extract the oldest in the target and the most recent date in the
	% reference
	if ~isempty(target_info)
		target_date = min([target_info.datenum]);
	else
		target_date = [];
	end
	
	if ~isempty(reference_info)
		reference_date = max([reference_info.datenum]);
	else
		reference_date = [];
	end
	
	
	% Handle the case where one does not exist
	if isempty(reference_date)
		older = 0;
	elseif isempty(target_date)
		older = 1;
	else % Compare the dates
		older = target_date < reference_date;
	end

end
