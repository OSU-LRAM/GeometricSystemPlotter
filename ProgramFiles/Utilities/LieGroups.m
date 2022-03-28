% LieGroups.m
% defines supported Lie Groups acting as configuration spaces
% this is a half-measure. groups should be classes; one day, update
% sysplotter to reflect the below capability.

classdef LieGroups
    % supported groups
    enumeration
        SE2, SO3
    end
    % group-specific helpers
    methods
        % dimensionality
        function dim = n_dim(group)
            % extend here for other groups
            if group == LieGroups.SE2 || group==LieGroups.SO3
                dim = 3;
            end
        end
        % matrix mapping
        function fn = mat_fn(group)
            if group == LieGroups.SE2
                fn = @vec_to_mat_SE2_lie;
            elseif group == LieGroups.SO3
                fn = @vec_to_mat_SO3;    
            end
        end
        % vector mapping
        function fn = vec_fn(group)
            if group == LieGroups.SE2
                fn = @mat_to_vec_SE2_lie;
            elseif group == LieGroups.SO3
                fn = @mat_to_vec_SO3;
            end
        end
        % optimizer
        function fn = optimizer_fn(group)
            if group == LieGroups.SE2
                fn = @optimize_coordinate_choice;
            elseif group == LieGroups.SO3
                fn = @optimize_so3_sys;
            end
        end
    end
end