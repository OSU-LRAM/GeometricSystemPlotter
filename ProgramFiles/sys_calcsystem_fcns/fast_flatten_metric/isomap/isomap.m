function [final_x, final_y, rv, D,EI] = isomap(x, y, springs, neutral_lengths, w, plot3d)
% Runs underlying MDS to finish Isomap
% x, y, springs, and neutral_lengths are from springs computation
% plot3d tells whether or not to plot 3d manifold
% w is weighting factor
	
    % Use spring lengths to set up sparse distance matrix
    n = numel(x);
    [D, G] = graph_matrix(springs, neutral_lengths, n);
    
%     Call Isomap 
    options.dims = 1:10;
    options.display = 0;
    options.verbose = 1;
    options.G = G;
    [Y, R, E] = isomap_fast(D, 'k', 8, options);

    % Call MDS
    weights = G~=0;
    weights = w*(1+weights)-(w-1);
%     [Y, R] = mdimscale(D, 1:3, 'stress', ones(size(weights)));
    
    % 3D Mesh plot
    if plot3d == 1
        a = Y.coords{3,1}(1,:);
        b = Y.coords{3,1}(2,:);
        c = Y.coords{3,1}(3,:);
        dx = .05;
        dy = .05;

        x_edge = [floor(min(a)):dx:ceil(max(a))];
        y_edge = [floor(min(b)):dy:ceil(max(b))];
        [A,B] = ndgrid(x_edge,y_edge);
        C = griddata(a,b,c,A,B);

%         figure;
% %         mesh(A,B,C);
%         surf(A,B,C);
%         axis equal;
    end
    
    new_coords = Y.coords{2,1};
    final_x = reshape(new_coords(1,:), size(x));
    final_y = reshape(new_coords(2,:), size(y));
    rv = R(2);
    
    % export extra info
    EI.A = A;
    EI.B = B;
    EI.C = C;
    
end