res = mopen ('E:\SciLabProg\ЧМ\Newton.csv', 'wt');

eps = 10^-9; i = 0; 
x1 = 0; x2 = 0;

y1 = 1; y2 = 0;
f1 = x1*x1 - x2 + 1; f2 = x1 - cos(%pi*x2/2);

mfprintf (res, 'I, X1, X2, F1, F2\n');
mfprintf (res, '%d, %f, %f, %f, %f\n', i, y1, y2, f1, f2);

while ( sqrt( (y1 - x1)^2 + (y2 - x2)^2 ) >= eps )
    x1 = y1; x2 = y2;
    
    i = i + 1;
    
    h = ( x2 - 1 + x1^2 - 2*x1*cos(%pi*x2/2) )/( -%pi*x1*sin(%pi*x2/2) - 1 );
    g = - x1 + cos(%pi*x2/2) - %pi*sin(%pi*x2/2)*h/2;
    
    y1 = x1 + g;
    y2 = x2 + h;
    
    f1 = y1^2 - y2 + 1;
    f2 = y1 - cos(%pi*y2/2);
    mfprintf (res, '%d, %f, %f, %f, %f\n', i, y1, y2, f1, f2);
end

disp ('x1:', x1); disp ('x2:', x2); 
disp ('f1:', f1); disp ('f2:', f2); 

mclose (res);
