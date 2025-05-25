function [Q, R] = qrgramschmidt(A)
A = [1 2 3 4; 1 3 5 7; 1 4 7 9; 1 5 8 11];
[m,n] = size(A);
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);
for k = 2:n
    R(1:k-1,k) = Q(1:m,1:k-1)'*A(1:m,k);
    z = A(1:m,k) - Q(1:m,1:k-1)*R(1:k-1,k);
    R(k,k) = norm(z);
    Q(1:m,k) = z/R(k,k);
end
Y = [30; 50; 66; 89];
beta = R\(Q'*Y);

disp('A=');
disp(A);
% Display Q and R as fractions
disp('Q=');
disp(rats(Q)); % Display Q as fractions

disp('R=');
disp(rats(R)); % Display R as fractions

disp('Vector beta:');
disp(rats(beta)); % Display beta as fractions
end
