function M_alpha = mass_matrix(geometry,physics,jointangles)
[A, ~,~,~,~, M_full, ~] = Inertial_connection_discrete(geometry,physics,jointangles);
M_alpha = mass_pull_back(M_full,A);
end