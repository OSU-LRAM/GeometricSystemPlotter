function P = interp_posdef(grid,posdef_values,evalpoint)


    % Make an eigen-decomposition of the positive-definite matrices
    [V, D] = cellfun(@(P) eig(P),celltensorconvert(posdef_values),'UniformOutput',false);
    
    % Interpolate the scale values
    D_i = cellfun(@(d) interpn(grid{:},d,evalpoint{:}),celltensorconvert(D));%,'UniformOutput',false);
    
    %%%% Interpolate the eigenvector rotation matrices
    
    
    %-----------------%
    
    % Convert the eigenvectors into a rotation matrix corresponding to minimally-rotated *version* This is a
    % hack at first, should incorporate ideas from paper below.
    % https://arxiv.org/pdf/1406.3361v1.pdf
    
    
    V = cellfun(@(v) flip_ith_column(v,1),V,'UniformOutput',false);
    
    
    % Corrected version should try all flips, take the logm (as below) for
    % each of them, and select the flip for each matrix that produces the
    % minimal rotation angle (will need to be careful then for sign flips,
    % etc)
    
    %----------------%
    
    % then take the matrix log of each of the rotation matrices
    Vcirc = cellfun(@(v) logm(v),V,'UniformOutput',false);
    
    % interpolate in the log space
    Vcirc_i = cellfun(@(v) interpn(grid{:},v,evalpoint{:}),celltensorconvert(Vcirc));%,'UniformOutput',false);
    
    % exponentiate the averaged matrix log
    V_i = expm(Vcirc_i); %cellfun(@(v) expm(v),Vcirc_i);%,'UniformOutput',false);
    
    
    %%%%%
    % construct the interpolated positive-definite matrix
    
    P = V_i * D_i;



end

function V = flip_ith_column(V,idx)

    V(:,idx) = -V(:,idx);

end