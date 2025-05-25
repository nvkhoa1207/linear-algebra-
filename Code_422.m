%4.2.2
% Define the matrix A
A = [2, 3, 4;
     4, 2, 5;
     5, 7, 2];

% Set tolerance and maximum iterations
tolerance = 1e-6; 
max_iter = 5; % Set maximum iterations to 5 

% Initialize variables
n = size(A, 1);% Dimension of the matrix
Ak = A; % Copy of A for iterative process [
eigenvalues = zeros(n, 1); % Preallocate space for eigenvalues 

% QR iteration to approximate eigenvalues
for k = 1:max_iter
    % Perform QR decomposition
    [Q, R] = qr(Ak);
    % Update Ak
    Ak = R * Q;

    % Check for convergence
    if norm(diag(Ak) - eigenvalues) < tolerance
        break;
    end

    % Update eigenvalues with diagonal elements
    eigenvalues = diag(Ak);
end

% Display the results
disp(diag(Ak));
