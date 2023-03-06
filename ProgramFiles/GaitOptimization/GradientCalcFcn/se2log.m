function xi = se2log(expXi)
    xi = zeros(3,1);
    % check if expXi is matrix or vector form
    if numel(expXi) == 3
        z_xy = expXi(1:2);
        z_theta = expXi(3);
    elseif numel(expXi) == 9
        z_xy = expXi(1:2,3);
        z_theta = atan2(expXi(2,1),expXi(1,1));
    end
    
    if abs(z_theta) < 1e-3
        xi(1:2) = z_xy;
        xi(3) = z_theta;
    else
        xi(1:2) = z_theta*....
            [sin(z_theta)/(2*(1-cos(z_theta))), 1/2;...
            -1/2, sin(z_theta)/(2*(1-cos(z_theta)))]*z_xy;
    
        xi(3) = z_theta;
    end
end