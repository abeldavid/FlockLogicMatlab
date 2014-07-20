function H2 = H2_norm(L, Q)

% Inputs: L (Laplacian matrix), Q (projection matrix
% onto n-1 dimensional space perpendicular to ones(n,1))

Lbar = Q*L*Q';
L0 = real(reduced_lyap(Lbar));
H2 = sqrt(trace(L0));