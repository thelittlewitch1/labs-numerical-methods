clear;

/*Условия задачи*/
A= [7.77, 0.27, - 0.29;
      1.15, - 6.22, 1.77;
      1.05, 4.52, 9.544];
b= [1.450; 1.050; - 1.310];
eps = 0.001;

/*Вспомоогательные фунции*/

/*Необходимые переменные*/
n = sqrt(length (A))
X = ones(1,n)';
p = zeros(1,n)';
counter = 0;
D = zeros(n,n); E = eye(n,n); invD = zeros(n,n);

/*Собственно тело программы*/
for i=1:1:n do
        D(i,i) = A(i,i);
        invD(i,i) = 1/D(i,i);
end
B = E - invD*A;
g = invD*b;

if ~(norm(B)<1) then
    disp ("I cannot solve this matrix equation, because the sufficient convergence condition is not met. Try another matrix.");
else
    while ~(norm(X - p) <= eps*(1-norm(B))/norm(B)) do
        p = X;
        X = B*p + g;
        counter = counter +1;
    end
    disp ("X:", X)
    disp ("counter:", counter)
    disp ("A*X:", A*X)
    disp ("b:", b)
    disp ("eps", eps)
    disp ("abs(A*X - b)", abs(A*X - b))
end;
