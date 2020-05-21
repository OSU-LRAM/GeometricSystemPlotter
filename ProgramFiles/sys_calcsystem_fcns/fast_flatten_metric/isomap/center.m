function cent = center(D, N, nl)
%CENTER Centers a matrix
% M is the matrix to be centered
% N is the original number of rows
% nl is the desired number of columns

%     cent = -.5*(D.^2 - sum(D'.^2)'*ones(1,nl)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
    cent =-0.5*((eye(N)-1/N*ones(N,N))*D.^2*(eye(N)-1/N*ones(N,N)));
end

