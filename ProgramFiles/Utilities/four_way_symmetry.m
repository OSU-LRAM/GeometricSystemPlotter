function [A,X1,X2,Y1,Y2,T1,T2,all_syms] = four_way_symmetry(alpha1, alpha2, Axa1, Axa2, Aya1, Aya2, At1,At2, a1, a2, lengthscale,mode)
% Take three connection vector fields for a three-link system, apply their necessary symmetry
% transformations, and average them. alpha1 and alpha2 are the joint angles where
% the fields are evaluated, a1 and a2 are the output locations
%
% Input fields should be the local gradient of position with respect to
% shape; the last section of this code negates this gradient to form the
% local connection

if ~exist('mode','var')
    mode = '3link';
end

% Make a list of fields
fields = {'x','y','theta'};

% Build up the lists of transformations
transforms = {'original','reflect_even','reflect_odd','rotate_180'};

%%%%%%%%%
% Define the transforms

switch mode
    
    case '3link'
        % Transforms for a three-link system (axes of symmetry are at 45 degrees to
        % shape parameters

        % Shape reflections
        f.a1.original = a1;
        f.a1.reflect_even = a2;
        f.a1.reflect_odd = -a2;
        f.a1.rotate_180 = -a1;

        f.a2.original = a2;
        f.a2.reflect_even = a1;
        f.a2.reflect_odd = -a1;
        f.a2.rotate_180 = -a2;

        % X field
        f.x1.original = Axa1;
        f.x1.reflect_even = -Axa2;
        f.x1.reflect_odd = Axa2;
        f.x1.rotate_180 = -Axa1;

        f.x2.original = Axa2;
        f.x2.reflect_even = -Axa1;
        f.x2.reflect_odd = Axa1;
        f.x2.rotate_180 = -Axa2;

        % Y field
        f.y1.original = Aya1;
        f.y1.reflect_even = Aya2;
        f.y1.reflect_odd = Aya2;
        f.y1.rotate_180 = Aya1;

        f.y2.original = Aya2;
        f.y2.reflect_even = Aya1;
        f.y2.reflect_odd = Aya1;
        f.y2.rotate_180 = Aya2;

        % Theta field
        f.theta1.original = At1;
        f.theta1.reflect_even = -At2;
        f.theta1.reflect_odd = -At2;
        f.theta1.rotate_180 = At1;

        f.theta2.original = At2;
        f.theta2.reflect_even = -At1;
        f.theta2.reflect_odd = -At1;
        f.theta2.rotate_180 = At2;
        
    case 'EvenOdd'
        
        % Shape reflections
        f.a1.original = a1;
        f.a1.reflect_even = a1;
        f.a1.reflect_odd = -a1;
        f.a1.rotate_180 = -a1;

        f.a2.original = a2;
        f.a2.reflect_even = -a2;
        f.a2.reflect_odd = a2;
        f.a2.rotate_180 = -a2;

        % X field
        f.x1.original = Axa1;
        f.x1.reflect_even = -Axa1;
        f.x1.reflect_odd = Axa1;
        f.x1.rotate_180 = -Axa1;

        f.x2.original = Axa2;
        f.x2.reflect_even = Axa2;
        f.x2.reflect_odd = -Axa2;
        f.x2.rotate_180 = -Axa2;

        % Y field
        f.y1.original = Aya1;
        f.y1.reflect_even = Aya1;
        f.y1.reflect_odd = Aya1;
        f.y1.rotate_180 = Aya1;

        f.y2.original = Aya2;
        f.y2.reflect_even = -Aya2;
        f.y2.reflect_odd = -Aya2;
        f.y2.rotate_180 = Aya2;

        % Theta field
        f.theta1.original = At1;
        f.theta1.reflect_even = -At1;
        f.theta1.reflect_odd = -At1;
        f.theta1.rotate_180 = At1;

        f.theta2.original = At2;
        f.theta2.reflect_even = At2;
        f.theta2.reflect_odd = At2;
        f.theta2.rotate_180 = At2;   
        
 case 'None'
        
        % Shape reflections
        f.a1.original = a1;
        f.a1.reflect_even = a1;
        f.a1.reflect_odd = a1;
        f.a1.rotate_180 = a1;

        f.a2.original = a2;
        f.a2.reflect_even = a2;
        f.a2.reflect_odd = a2;
        f.a2.rotate_180 = a2;

        % X field
        f.x1.original = Axa1;
        f.x1.reflect_even = Axa1;
        f.x1.reflect_odd = Axa1;
        f.x1.rotate_180 = Axa1;

        f.x2.original = Axa2;
        f.x2.reflect_even = Axa2;
        f.x2.reflect_odd = Axa2;
        f.x2.rotate_180 = Axa2;

        % Y field
        f.y1.original = Aya1;
        f.y1.reflect_even = Aya1;
        f.y1.reflect_odd = Aya1;
        f.y1.rotate_180 = Aya1;

        f.y2.original = Aya2;
        f.y2.reflect_even = Aya2;
        f.y2.reflect_odd = Aya2;
        f.y2.rotate_180 = Aya2;

        % Theta field
        f.theta1.original = At1;
        f.theta1.reflect_even = At1;
        f.theta1.reflect_odd = At1;
        f.theta1.rotate_180 = At1;

        f.theta2.original = At2;
        f.theta2.reflect_even = At2;
        f.theta2.reflect_odd = At2;
        f.theta2.rotate_180 = At2;  
        
    otherwise
        
        error('Invalid symmetry averaging mode')
        
end

% Interpolate the 4 transformed fields for each component at the joint
% angles
for i = 1:length(fields)
			
	for j = 1:2 % column-number
	
		for k = 1:length(transforms)
			
			all_syms.(fields{i}){k,j} = interpn(...
				alpha1 ...f.alpha1.(transforms{k})...
				,alpha2 ...,f.alpha2.(transforms{k})...
				,f.([fields{i} num2str(j)]).(transforms{k})...
				,f.a1.(transforms{k})...
				,f.a2.(transforms{k})...,a1,a2 ...
				);
		

		
		end
		
		% Average the fields
		dimfields = numel(size((all_syms.(fields{i}){1,1})));
		averaged_fields.(fields{i}){j} = sum(cat(dimfields+1,all_syms.(fields{i}){:,j}),dimfields+1)/size(all_syms.(fields{i}),1);
		
	end
	
end

% Extract the output
A = {-lengthscale*averaged_fields.x{1} -lengthscale*averaged_fields.x{2};
	-lengthscale*averaged_fields.y{1} -lengthscale*averaged_fields.y{2};
	-averaged_fields.theta{1} -averaged_fields.theta{2}};

% Extract the output
X1 = -lengthscale*averaged_fields.x{1};
X2 = -lengthscale*averaged_fields.x{2};

Y1 = -lengthscale*averaged_fields.y{1};
Y2 = -lengthscale*averaged_fields.y{2};

T1 = -averaged_fields.theta{1};
T2 = -averaged_fields.theta{2};
end