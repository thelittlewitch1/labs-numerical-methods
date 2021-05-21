function res = isSA (M)
    flag = %T
    if (size(M, 1) == size(M, 2)) then
        n = size(M, 1)
        for i = 1:n
            for j = 1:n 
                if ( M(i,j) ~= conj(M(j,i)) )
                    flag = %F;
                end;
            end
        end
    else
        flag = %F;
    end
    res = flag;
endfunction

function [u0, lambd] = fun (A, u, eps)
    lambd0 = 0;       
    lambd1 = max ( A  * u );

    k = 0;
    //disp('---------')
    //disp ('lambda = ', lambd1)
    //disp ( '( A^(k+1) ) * u = ', ( A^(k+1) )*u )
    while (abs(lambd0 - lambd1) > eps)
        lambd0 = lambd1;
        k = k+1;
        lambd1 = max(  ( A^(k+1) ) * u  )/max( A^k * u)
        //disp ('lambda = ', lambd1)
        //disp ( '( A^(k+1) ) * u = ', ( A^(k+1) )*u )
    end
    //disp('---------')
    u0 = ( A^(k+1) )*u;
    lambd = lambd1;
endfunction

A = [4, 1, 1;
        1, 2, 1;
        1, 1, 3;]
u = [1;0; 1]
eps = 10^-3

if ( isSA(A) && (size(A,2) == size(u, 1) ) ) then
    disp("A:", A, "u:", u);
    [u0, lambd] = fun (A, u, eps);
    disp("total max lambda:", lambd, 'corresponding vec: ', u0);

    disp("A^-1:", A^-1, "u:", u);
    [u0, lambd] = fun (A^-1, u, eps)
    disp("total min lambda:", lambd, 'corresponding vec: ', u0);
else
    disp ('Matrix A: ', A, 'is not self-adjoint (Hermitian) matrix', 'or something else is wrong.')
end

