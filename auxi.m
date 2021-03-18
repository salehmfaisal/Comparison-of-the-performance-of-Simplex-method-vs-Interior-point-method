
 function Ap = auxi(A,b)
 [m,~] = size(A);
 for i = 1:m
     if b(i) == 0
         b(i) = 1;
     else
        b(i) = sign(b(i));
     end
 end
 b1 = diag(b);
 Ap = [A b1];
 
 