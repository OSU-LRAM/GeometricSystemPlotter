function [g_end_orig,g_end_opt, cost_end,BVI_orig, cBVI_orig, BVI_opt,cBVI_opt] = extract_displacement_and_cost(datafile)
% Extract the displacement and cost data from a sysf_...shchf_....mat file

% Load the target file
load(datafile,'p')

% Prime arrays to hold the net displacement (in original and optimal
% coordinates) and cost from each shape change in the file. p.G_locus_full is
% single-level cell array of structures, each of which holds the
% information for one gait (with all segments concatenated)
g_end_orig = zeros(numel(p.G_locus_full),3);
g_end_opt = g_end_orig;
cost_end = zeros(numel(p.G_locus_full,1)); % If distance metric was not specified, euclidean metric in the parameters was assumed
BVI_orig = g_end_orig;
cBVI_orig = g_end_orig;
BVI_opt = g_end_orig;
cBVI_opt = g_end_orig;

% Loop over each shape change
for i = 1:numel(p.G_locus_full)
	
	% Extract the final values for the relevant parameters
	g_end_orig(i,:) = p.G_locus_full{i}.G(end,:); 
	g_end_opt(i,:) = p.G_locus_full{i}.G_opt(end,:); 
	cost_end(i) = p.G_locus_full{i}.S(end);
    
    BVI_orig(i,:) = p.G_locus_full{i}.bvi(end,:);
    BVI_opt(i,:) = p.G_locus_full{i}.bvi_opt(end,:);
    cBVI_orig(i,:) = p.cBVI{i};
    cBVI_opt(i,:) = p.cBVI_opt{i};
    
end