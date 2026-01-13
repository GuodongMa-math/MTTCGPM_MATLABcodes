% Matlab Model by Jianghua Yin (Apr.,2019,Nanning)
% Copyright (C) 2019 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem8(x)

n = length(x);
out = zeros(n,1);
out(1) =4* x(1)*(x(1)^2+x(2)^2)-4;
for i=2:n-1
    out(i) =4*x(i)*(x(i-1)^2+x(i)^2)+4*x(i)*(x(i+1)^2+x(i)^2)-4;
end
out(n) =4* x(n)*(x(n-1)^2+x(n)^2);


    
        
    
    


