% Matlab Model by Jianghua Yin (Nov.,2017,Laibin)
% Copyright (C) 2017 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem7(x)

n = length(x);
out = zeros(n,1);
for i=1:n
    out(i)=(i/n)*exp(x(i))-1;
end


