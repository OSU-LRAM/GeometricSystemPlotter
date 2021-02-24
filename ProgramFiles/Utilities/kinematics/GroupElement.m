% generalization of a Lie Group element in SE(2) or SE(3)
classdef GroupElement
    % relevant representations
    properties
        is_planar;
        vector;
        matrix;
    end
    
    % static transforms, for convenience
    properties (Constant)
        % planar
        LINEAR = GroupElement([1 0 0]);
        ROTARY = GroupElement([0 0 1]);
        % 3d
        LINEAR_X = GroupElement([1 0 0 0 0 0]);
        LINEAR_Y = GroupElement([0 1 0 0 0 0]);
        LINEAR_Z = GroupElement([0 0 1 0 0 0]);
        ROTARY_X = GroupElement([0 0 0 1 0 0]);
        ROTARY_Y = GroupElement([0 0 0 0 1 0]);
        ROTARY_Z = GroupElement([0 0 0 0 0 1]);
    end
    
    methods
        % produces a GroupElement with vector in 2D or 3D
        function obj = GroupElement(rep)
            % sanitize inputs
            if nargin == 0
                rep = [0 0 0];
            end
            if isrow(rep)
                rep = rep';
            end
            % handle matrix construction
            if size(rep,2) > 1
                % rep is already a matrix; not setting vector
                obj.matrix = rep;
                warning("GroupElement constructed from matrix; no vector representation available")
                % do planar detection
                if length(obj.matrix) == 3
                    obj.is_planar = true;
                elseif length(obj.matrix) == 4
                    obj.is_planar = false;
                else
                    error('GroupElement must be in SE(2) or SE(3)')
                end
                return
            else
                % construct matrix from vector
                % do planar detection
                if length(rep) == 3
                    obj.is_planar = true;
                elseif length(rep) == 6
                    obj.is_planar = false;
                else
                    error('GroupElement must be in SE(2) or SE(3)')
                end
                if obj.is_planar
                    % 2d case is trivial
                    obj.matrix = [cos(rep(3)) -sin(rep(3)) rep(1);...
                                  sin(rep(3)) cos(rep(3)) rep(2);...
                                  0 0 1];
                else
                    % 3d case needs more complexity
                    obj.matrix = [GroupElement.rot_mat_3d(rep(4:end)), rep(1:3); [0 0 0 1]];
                end
                % save vector for easy access
                obj.vector = rep;
            end
        end
        
        % log to convert to Algebra element
        function lae = to_AlgebraElement(obj)
            lae = AlgebraElement(logm(obj.matrix));
        end
        
        % group actions (convenience functions)
        function obj = rightAction(obj, other)
            obj.matrix = obj.matrix * other.matrix;
            obj.vector = [];
            warning("GroupElement constructed from matrix; no vector representation available")
        end
        function obj = leftAction(obj, other)
            obj.matrix = other.matrix * obj.matrix;
            obj.vector = [];
            warning("GroupElement constructed from matrix; no vector representation available")
        end
    end
    
    methods(Static)
        % produces a 3D rotation matrix, in YPR order
        function R = rot_mat_3d(rot_vec)
            Rx = [1 0 0;
                  0 cos(rot_vec(1)) -sin(rot_vec(1));
                  0 sin(rot_vec(1)) cos(rot_vec(1))];
            Ry = [cos(rot_vec(2)) 0 sin(rot_vec(2));
                  0 1 0;
                  -sin(rot_vec(2)) 0 cos(rot_vec(2))];
            Rz = [cos(rot_vec(3)) -sin(rot_vec(3)) 0;
                  sin(rot_vec(3)) cos(rot_vec(3)) 0;
                  0 0 1];
            R = Rz * Ry * Rx;
        end
    end
end