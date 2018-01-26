function L_t = extract_length_over_time(stretch_const,L,backbone_type,backbone_section,shapechange)
% extract a time array describing the total length of the snake backbone
% over one cycle. 
%
% stretch_const ---- changeable value in lambda(a,s) function. usually between
%                    0 and 1. 
% L ---------------- Total length of backbone (m). Should always be 1, but
%                    left as an input for future modifications. 
% backbone_type ---- For current version of sysplotter, backbone types are
%                    serpenoid and piecewise. Enter as a string. 
% backbone section - A percentage between 0 and 1 for which this function
%                    will evaluate the length starting at the head and
%                    traveling backwards. For example, 0.5 will evaluate
%                    between the head and midpoint of the backbone. 
% shapechange ------ Which shapechange data to investigate. Include shchf
%                    but no file extension. 

b = stretch_const;

if strcmp(backbone_type,'serpenoid')
    curvdef = str2func(['curv_serpenoid_stretch_b_',num2str(b*1000)]);
    zfile = ['sysf_serpenoid_extendable_',num2str(b*1000),'__',shapechange,'.mat']; % shchf_circle_6_dphi
elseif strcmp(backbone_type,'piecewise')
    curvdef = str2func(['curv_piecewise_const_stretch_b_',num2str(b*1000)]);
    zfile = ['sysf_piecewise_const_stretch_b_',num2str(b*1000),'__',shapechange,'.mat'];
%     curvdef = str2func('curv_piecewise_const_stretch_abs_b_100');
%     zfile = ['sysf_piecewise_const_stretch_abs_b_100__',shapechange,'.mat']; % override for testing
else
    error('Unknown backbone type.')
end


% load target file
load(zfile,'p');

% extract time array and alpha inputs
t = p.time{1}{1};
a = p.phi_fun{1}{1}(t);

% for this particular system, go through each time step
% length_end = floor(backbone_section*numel(t));
L_t = zeros(numel(t),1);

% for each time step in the gait
for j = 1:numel(t)
    
    cparams = a(j,:)';
    
    % extract backbone
    [h,~] = backbone_from_stretchable_curvature(curvdef,cparams,L);
    
    % define x, y, th for this time
    s_loc = h(linspace(L*-.5,L*.5,numel(t)))';
    
    % extract total length by adding each step along backbone in specified range
    for k = floor(size(s_loc,1)-backbone_section*size(s_loc,1))+2:size(s_loc,1)
        L_t(j) = L_t(j) + sqrt((s_loc(k,1)-s_loc(k-1,1))^2 + (s_loc(k,2)-s_loc(k-1,2))^2);
    end
    
end


end