function [D, G] = graph_matrix(springs, neutral_lengths, n)
%GRAPH_MATRIX Computes the matrix of distances given springs and lengths
%   D: all distances, inferred using Dijkstra
%   G: only neighbor distances from spring lengths
    G = inf*ones(n);
    for i=1:n
        G(i,i)=0;
    end
    idx = sub2ind(size(G), springs(:,1), springs(:,2));
    G(idx) = neutral_lengths;
    idx = sub2ind(size(G), springs(:,2), springs(:,1));
    G(idx) = neutral_lengths;
    D = dijkstra(G,1:n);
end

