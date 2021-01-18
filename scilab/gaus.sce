clear;

/*Условия задачи*/
/*A= [7.77, 0.27, - 0.29;
      1.15, - 6.22, 1.77;
      1.05, 4.52, 9.544];
b= [1.450; 1.050; - 1.310];*/

A=read('E:\SciLabProg\A.txt', 3,3);
b=read('E:\SciLabProg\b.txt', 3,1);


/*Вспомоогательные фунции*/
function [z] = Exchange(C,i)
    k=i+1;
    
    while C(k,i)==0
        k=k+1;
    end;
    for j=1:size(C,1)
        s=C(i,j);
        C(i,j)=C(k,j);
        C(k,j)=s;
    end;
    z=C;
endfunction

function [z] = Simplex(A,b)
    N=size(A,1); //определение числа уравнений системы
    C=cat(2,A,b); //создание расширенной матрицы системы
    for i=1:N-1
        if C(i,i)==0
            C=Exchange(C,i);
        end;
        for j=0:N
            C(i,N+1-j)=C(i,N+1-j)/C(i,i);
        end;
        for m=i+1:N
            alpha=C(m,i);
            for j=i:N+1
                C(m,j)=C(m,j)-alpha*C(i,j);
            end;
        end;
    end;
    C(N,N+1)=C(N,N+1)/C(N,N);
    C(N,N)=1;
    z=C;
endfunction

function [z] = Gauss(A,b)
    C=Simplex(A,b);
    N=size(A,1);
    v(N)=C(N,N+1);
    for j=1:N-1
        s=0;
        for k=0:j-1
            s=s+C(N-j,N-k)*v(N-k);
        end;
        v(N-j)=(C(N-j,N+1)-s)/C(N-j,N-j);
    end;
    z=v';
endfunction

/*Собственно тело программы*/
x = Gauss(A, b)'; 

disp ("A:", A)
disp ("b:", b)
disp ("X:", x)
disp ("A*X:", A*x)

