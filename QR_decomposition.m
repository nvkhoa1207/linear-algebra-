function [Q, R] = GramSchmidtQR(A)
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
A_example = [0 4 3; -1 2 4; 1 0 0]; 
[Q, R] = GramSchmidtQR(A_example);

disp('Q matrix from GramSchmidtQR:');
disp(Q);
disp('R matrix from GramSchmidtQR:');
disp(R);