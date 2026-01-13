% Matlab Model by Jianghua Yin (Nov.,2017,Laibin)
% Copyright (C) 2017 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem6(x)

n = length(x);
out = zeros(n,1);
out(1) = x(1)-exp(cos((x(1)+x(2))/2));
for i=2:n-1
    out(i) = x(i)-exp(cos((x(i-1)+x(i)+x(i+1))/(n+1)));
end
out(n) = x(n)-exp(cos((x(n-1)+x(n))/(n+1)));


    
        
    
    


