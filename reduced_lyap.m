function X = reduced_lyap(Lbar)
%REDUCED_LYAP function to solve the Lyapunov equation 
%Lbar.X + X.Lbar* = I for reduced Laplacian matrices Lbar.
%   X = REDUCED_LYAP(Lbar) returns the symmetric matrix X which solves the
%   Lyapunov equation for the given reduced Laplacian matrix Lbar.
%
%   This routine first converts Lbar to Schur form, since reduced Laplacian
%   matrices are equivalent up to unitary transforms. Then, the Lyapunov
%   equation is solved for the upper triangular Schur matrix. 

[N1 N2] = size(Lbar);

if (N1 ~= N2)
    disp('ERROR: Input must be a square matrix')
else        
    T = schur(Lbar, 'complex');

    X = zeros(N1,N1);

    for r = N1:-1:1
        for c = N1:-1:r

            if r == c
                rhs = 1;
            else
                rhs = 0;
            end

            for i = (r+1):N1
                rhs = rhs - T(r,i)*X(i,c);
            end

            for i = (c+1):N1
                rhs = rhs - X(r,i)*conj(T(c,i));
            end

            X(r,c) = rhs/(T(r,r) + conj(T(c,c)));
            X(c,r) = conj(X(r,c));
            
        end
    end
end