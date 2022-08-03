function ms = MeanSquare2Norm( array1, array2 )
%MEANSQUARE2NORM Summary of this function goes here
%   Detailed explanation goes here


ms = (array1 - array2).^2;
ms = sum(ms(:)) / numel(array1);



end

