%Classical GramSchmidt
function [Q, R] = group10_topic6(A)
[m, n] = size(A);

if m < n
    error('Number of rows must be greater than or equal to number of columns.');
end

Q = zeros(m, n);
R = zeros(n, n);

for j = 1:n
    v_j = A(:, j); 
    for i = 1:j-1
        R(i, j) = Q(:, i)' * A(:, j); 
        v_j = v_j - R(i, j) * Q(:, i);
    end
    R(j, j) = norm(v_j);
    if R(j, j) < eps 
        error('Matrix does not have full column rank or is ill-conditioned.');
    end
    Q(:, j) = v_j / R(j, j);
end
end

%Modified GramSchmidt
function [Q, R] = modifiedGramSchmidtQR(A)

[m, n] = size(A);

if m < n
    error('Number of rows must be greater than or equal to number of columns.');
end

Q = zeros(m, n);
R = zeros(n, n);
V = A;

for i = 1:n
    R(i, i) = norm(V(:, i));
    if R(i,i) < eps
        error('Matrix does not have full column rank or is ill-conditioned.');
    end
    Q(:, i) = V(:, i) / R(i, i);
    for j = i+1:n
        R(i, j) = Q(:, i)' * V(:, j);
        V(:, j) = V(:, j) - R(i, j) * Q(:, i);
    end
end
end

%4.1.1
function beta = calculateLinearRegressionQR(A_input, Y_input)

[m,n] = size(A_input);

% Initialize R and Q matrices
R = zeros(n, n);
Q = zeros(m, n);

% Gram-Schmidt QR Decomposition
R(1,1) = norm(A_input(:,1));
Q(:,1) = A_input(:,1) / R(1,1);

for k = 2:n
    R(1:k-1,k) = Q(:,1:k-1)' * A_input(:,k);
    z = A_input(:,k) - Q(:,1:k-1) * R(1:k-1,k);
    R(k,k) = norm(z);
    Q(:,k) = z / R(k,k);
end

% Solve for beta using R*beta = Q'*Y
beta = R \ (Q' * Y_input);

disp('A=');
disp(A_input);

disp('Q=');
disp(rats(Q)); % Display Q as fractions

disp('R=');
disp(rats(R)); % Display R as fractions

end

%4.2.2
function eigenvalues = calculateEigenvaluesQR(A, tolerance, max_iter)
% Initialize variables
n = size(A, 1);
Ak = A;
eigenvalues = zeros(n, 1);

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
    % Update eigenvalues with diagonal elements for the next iteration's comparison
    eigenvalues = diag(Ak);
end

% The final approximated eigenvalues are the diagonal elements of Ak
eigenvalues = diag(Ak);
end



A_example = [0 4 3; -1 2 4; 1 0 0]; 
[Q, R] = GramSchmidtQR(A_example);

%CGS
disp('Q matrix from GramSchmidtQR:');
disp(Q);
disp('R matrix from GramSchmidtQR:');
disp(R);

%MGS
[Q1, R1] = modifiedGramSchmidtQR(A_example);
disp('Q matrix from modifiedGramSchmidtQR:');
disp(Q1);
disp('R matrix from modifiedGramSchmidtQR:');
disp(R1);

%4.1.1
disp('4.1.1');
A_411 = [1 2 3 4; 1 3 5 7; 1 4 7 9; 1 5 8 11];
Y = [30; 50; 66; 89];
beta = calculateLinearRegressionQR(A_411, Y);

disp('Vector beta:');
disp(rats(beta)); % Display beta as fractions

%4.2.2
disp('4.2.2');
% Define the matrix A
A_422 = [2, 3, 4; 4, 2, 5; 5, 7, 2];
tol = 1e-6;
maxIter = 5;
ans = calculateEigenvaluesQR(A_422, tol, maxIter);
disp('Approximated Eigenvalues:');
disp(ans);
