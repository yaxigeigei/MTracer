function ms = XMeanSquare( fixed, moving, msLength, focus )
%MEANSQUARE2 Summary of this function goes here
%   Detailed explanation goes here
    
    
    
    
    fLength = length(fixed);
    mLength = length(moving);
    fullLength = fLength + mLength - 1;
    
    if nargin < 4
        msLength = fullLength;
        focus = ceil(msLength/2);
    end
    
    
    ms = zeros(1, msLength);
    iBegin = ceil(fullLength/2)-ceil(msLength/2)+1;
    range = (iBegin : iBegin+msLength-1) + focus-ceil(fLength/2); % cautious
    
    mskLength = fLength + 2 * (mLength - 1);
    mskF = cat(2, zeros(1, mLength-1), ones(1, fLength), zeros(1, mLength-1));
    
    
    for i = range
        
        mskM = zeros(1, mskLength);
        mskM(i : i-1+mLength) = 1;
        mskShare = mskM .* mskF;
        
        mskIndex = find(mskShare);
        mIndex = mskIndex - (i-1);
        fIndex = mskIndex - (mLength-1);
        
        ms(i-range(1)+1) = MeanSquare2Norm(fixed(fIndex), moving(mIndex));
        
    end
    



end





