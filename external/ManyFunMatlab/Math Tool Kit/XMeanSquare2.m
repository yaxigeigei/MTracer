function ms = XMeanSquare2(a1, a2, varargin)
%MEANSQUARE2 Summary of this function goes here
%   Detailed explanation goes here
    
    [ h1, w1 ] = size(a1);
    [ h2, w2 ] = size(a2);
    ms = zeros(h1+h2-1, w1+w2-1);
    
    p = inputParser();
    p.addParameter('Focus', [ ceil(h1/2), ceil(w1/2) ], @isnumeric);
    p.addParameter('OutputSize', [], @isnumeric);
    
    p.parse(varargin{:});
    focus = p.Results.Focus;
    outputSize = p.Results.OutputSize;
    
    if isempty(outputSize)
        outputSize = size(ms);
        rBeginInd = 1;
        rEndInd = outputSize(1) - 1;
        cBeginInd = 1;
        cEndInd = outputSize(2) - 1;
    else
        rBeginInd = ceil((h2-1)/2) + focus(1) - ceil(outputSize(1)/2);
        rEndInd = rBeginInd + outputSize(1) - 1;
        cBeginInd = ceil((w2-1)/2) + focus(2) - ceil(outputSize(2)/2);
        cEndInd = cBeginInd + outputSize(2) - 1;
    end
    
    [ rr1, rr2 ] = XIndexing(a1(:,1), a2(:,1));
    [ cr1, cr2 ] = XIndexing(a1(1,:), a2(1,:));
    
    a1 = double(a1);
    a2 = double(a2);
    
    for i = rBeginInd : rEndInd
        for j = cBeginInd : cEndInd
            ms(i,j) = MeanSquare2Norm(a1(rr1(i,1):rr1(i,2), cr1(j,1):cr1(j,2)), ...
                a2(rr2(i,1):rr2(i,2), cr2(j,1):cr2(j,2)));
        end
    end
    
    ms = ms(rBeginInd : rEndInd, cBeginInd : cEndInd);
    
    
end

