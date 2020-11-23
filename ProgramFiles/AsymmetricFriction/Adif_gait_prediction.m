function [endpoint] = Adif_gait_prediction(s, amplitude, offset)
%ADIF_GAIT_PREDICTION Summary of this function goes here
%   Detailed explanation goes here

A = get_A_of_alpha(s, grid).difference;

Ax_dif = A(:,1);
Ay_dif = A(:,2);
At_dif = A(:,3);

point_a = offset - amplitude;

point_a = a_grid(14)
point_b = a_grid(18)
da = a_grid(2) - a_grid(1)

dx = trapz(Ax_dif(14:18))*da
dy = trapz(Ay_dif(14:18))*da
dth = trapz(At_dif(14:18))*da

% compare to a gait that goes from point_a to point_b:
off = (point_a + point_b)/2;
amp = abs(point_a - point_b)/2;
gait = generate_1D_gait(amp,off,0);
sol = asym_solve_gait(s, gait, system_map);

displacement = deval(sol, 2*pi)'

g_cerc_vector = [dx dy dth];
g_circ = [0 -g_cerc_vector(3) g_cerc_vector(1); ...
          g_cerc_vector(3) 0 g_cerc_vector(2); ...
          0 0 0];
      % TODO: make this a function like vectomatSE2
      % vec_to_mat_se2 (v for vector)
      
Adif_prediction = mat_to_vec_SE2(expm(g_circ))

error = Adif_prediction - displacement
proport_error = error ./ displacement

end

