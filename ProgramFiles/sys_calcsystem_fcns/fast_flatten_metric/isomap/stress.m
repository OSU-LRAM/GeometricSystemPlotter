function S = stress(D, new_x, new_y)
%STRESS Computes the normalized stress

    D_new = pdist2([new_x(:) new_y(:)], [new_x(:) new_y(:)]);
    S = sqrt(sum(sum(triu(D-D_new).^2)) / sum(sum(triu(D).^2)));
end

