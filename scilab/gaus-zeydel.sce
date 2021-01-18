clear;

function res = iscollomzero (a,n)
    buf = 0; flag = %T
    for i=1:1:n do
        buf=0;
        for j=1:1:n
            buf = buf + a(i,j);
        end
        if buf == 0 then flag = %F end
    end
    res = flag;
endfunction

function res = isstrzero (a,n)
    buf = 0; flag = %T
    for i=1:1:n do
        buf=0;
        for j=1:1:n
            buf = buf + a(j,i);
        end
        if buf == 0 then flag = %F end
    end
    res = flag;
endfunction

/*Условия задачи*/
A= [7.77, 0.27, - 0.29;
      1.15, - 6.22, 1.77;
      1.05, 4.52, 9.544];
b= [1.450; 1.050; - 1.310];
eps = 0.001;


/*Вспомоогательные фунции*/
function [res] = isConvergent (a, n)
    flag = %T
    for i=1:1:n
        summ = 0;
        for j=1:1:n
            if j <> i then
                summ = summ + abs(a(i,j));
            end
        end
        if abs(a(i,i)) < summ
            then flag = %F;
        end
    end
    res = flag;
endfunction

function [res] = curdif (a, x, p, n, i)
    buf = 0;
    for j=1:1:i-1
        buf = buf - a(i,j)*x(j);
    end
    for j=i+1:1:n
        buf = buf - a(i,j)*p(j);
    end
    res = buf;
endfunction

function [res] =  isDiagContainZero (a, n)
    buf = 1; 
    for i=1:1:n do
        buf = buf * a(i,i);
    end
    res = (buf == 0);
endfunction

function [res] = curEps (x,p)
    m = abs(x(1) - p(1));
    for i=2:1:n
        m = max (m, abs(x(i)-p(i)));
    end
    res = m
endfunction

/*Необходимые переменные*/
n = sqrt(length (A))
X = ones(1,n)';
p = ones(1,n)'*0;
mprew = 0; 
counter = 0;

while isDiagContainZero(A,n) do
    for i=1:1:n
        if A(i,i) == 0 then
            if (i>1) then
                A = [A(:,1:i-2), A(:, i), A(:, i-1), A(:, i+1:n)]
            else
                A = [A(:,n), A(:, 1:n-1)];
            end
        end
    end
end

/*Собственно тело программы*/
if ~isConvergent (A, n) then
    disp ("I cannot solve this matrix equation, because there is no predominance of diagonal elements. Try another matrix.");
else
    while (curEps(X,p)>eps) do
        if mprew  == curEps(X,p) then
            disp ("I cannot solve this matrix equation.");
            abort;
        else
            mprew = curEps(X,p)
            p = X;

            for i=1:1:n do
                X(i) = ( b(i) +  curdif (A, X, p, n, i) ) / A(i,i)
            end
        end
        counter = counter +1;
    end
    disp ('A:', A)
    disp ('b:', b)
    disp ('eps', eps)
    disp ("counter:", counter)
    disp ("X:", X)
    disp ("A*X:", A*X)
    disp ("abs(A*X - b)", abs(A*X - b))
end;
