function [V,D] = jacobi(A,tol)
%JACOBI Compute the eigenvalues and eigenvectors using Jacobi
%   [V,D] = JACOBI(A,TOL) computes the eigenvectors V and eigenvalues D using the Jacobi
%   method. TOL is the tolerance used as a stopping criterion.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)

    V  = eye(size(A));
    while max(max(abs(triu(A, 1)))) > tol*max(abs(diag(A)))
        for i = 1:size(A, 1)
            for j = i+1:size(A, 2)
                theta = 0.5*atan(2*A(i,j)/(A(j,j)-A(i,i)));
                J = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                A([i,j],:) = J'*A([i,j],:);
                A(:,[i,j]) = A(:,[i,j])*J;
                V(:,[i,j]) = V(:,[i,j])*J;
            end
        end
    end
    D = diag(diag(A));
end
