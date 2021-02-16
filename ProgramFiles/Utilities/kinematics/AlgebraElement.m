% generalization of a Lie Algebra element in SE(2) or SE(3)
% almost a CC of GroupElement; copying and pasting code (bad practice)
classdef AlgebraElement
    properties
        is_planar;
        vector;
        matrix;
    end
    
    methods
        % construct LA element
        function obj = AlgebraElement(rep)
            % sanitize inputs
            if nargin == 0
                rep = [0 0 0];
            end
            if isrow(rep)
                rep = rep';
            end
            if size(rep,2) > 1
                % rep is already a matrix
                if length(rep) == 3
                    obj.is_planar = true;
                    obj.vector = [rep(1,3) rep(2,3) rep(2,1)]';
                elseif length(rep) == 4
                    obj.is_planar = false;
                    obj.vector = [rep(1,4) rep(2,4) rep(3,4) rep(3,2) rep(1,3) rep(2,1)];
                else
                    error('AlgebraElement must be in SE(2) or SE(3)')
                end
                obj.matrix = rep;
                return
            else
                % construct matrix
                % do planar detection
                if length(rep) == 3
                    obj.is_planar = true;
                elseif length(rep) == 6
                    obj.is_planar = false;
                else
                    error('AlgebraElement must be in SE(2) or SE(3)')
                end
                if obj.is_planar
                    obj.matrix = [0 -rep(3) rep(1);
                                  rep(3) 0 rep(2);
                                  0 0 0];
                else
                    obj.matrix = [0 -rep(6) rep(5) rep(1);
                                  rep(6) 0 -rep(4) rep(2);
                                  -rep(5) rep(4) 0 rep(3);
                                  0 0 0 0];
                end
            end
            obj.vector = rep;
        end
        % exponentiate to convert to group element
        function lge = to_GroupElement(obj, time)
            if ~exist('time','var')
                time = 1;
            end
            lge = GroupElement(expm(time * obj.matrix));
        end
    end
end