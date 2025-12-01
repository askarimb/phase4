pe = 6*p;
Lt = place(A', C', pe);   % Lt = L^T
L  = Lt.';                % 4x2 observer gain matrix

% 3) Verify that the observer poles are correctly placed
lambda_obs = eig(A - L*C);
