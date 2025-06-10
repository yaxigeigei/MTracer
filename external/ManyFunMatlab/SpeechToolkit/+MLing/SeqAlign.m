function [aa, steps, tscore] = SeqAlign(s1, s2, varargin)
% Globally align two sequences using Needleman-Wunsch algorithm
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/44617-sequence-alignment-with-arbitrary-steps
% 
%   [aa, steps, tscore] = SeqAlign(s1, s2)
%   [aa, steps, tscore] = SeqAlign(..., 'ScoreFunc', @scoreFunc)
%   [aa, steps, tscore] = SeqAlign(..., 'StepMatrix', [1 0; 0 1; 1 1])
% 
% Inputs
%   s1, s2              Two sequences to align.
%                       1) If s1 and s2 are vectors of char, they will be aligned at the level of 
%                          individual characters. Example: 'this is a sequence'.
%                       2) If s1 and s2 are vectors of string objects, they will be aligned at the 
%                          level of individual strings. Example: ["this" "is" "a" "sequence"].
%                       3) Aligning sequences of other datatypes is possible if elements can be 
%                          indexed with parentheses and support == operator, e.g. numeric sequences. 
%   'ScoreFunc'         A function handle that receives two inputs x and y, and returns a numeric 
%                       score. The default private function returns/scores 1 for match (x==y), 
%                       -1 for mismatch (x~=y) or gap (x or y is empty).
%   'StepMatrix'        A k-by-2 matrix that specify allowed steps. The default matrix shown above 
%                       implements the classical Needleman-Wunsch algorithm. 
% Outputs
%   aa                  Aligned sequences. This is a 2-by-n char or string array where the top row 
%                       is the aligned sequence of s1, bottom row s2. 
%                       For example:
%                       ['DID THOSE THIEVE-S TAKE THIRTY JEWELS';
%                        'DID THOSE J--EWELS TAKE T-IN-Y-------']
%                       Depending on the datatype, gaps are represented by '-' for char, "" for string, 
%                       NaN for numeric values, or duplicating the last element before the gap in case
%                       of other datatypes.
%   steps               Corresponding steps. This is a 2-by-n numeric array where the top row contain 
%                       steps to align s1, and the bottom for s2. Zeros represent gaps. The steps for 
%                       the above example would be:
%                       [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
%                        1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 0 0 0 0 0 0 0]
%   tscore              Total score of the alignment.
% 

% Options
p = inputParser();
p.addParameter('ScoreFunc', @scoreFunc, @(x) isa(x, 'function_handle'));
p.addParameter('StepMatrix', [1 0; 0 1; 1 1], @(x) size(x,2)==2);
p.parse(varargin{:});
fScore = p.Results.ScoreFunc;
S = p.Results.StepMatrix;

s1 = s1(:)';
s2 = s2(:)';
n = length(s1);
m = length(s2);

% First, retrieve the matrices V and T with best values and best steps, respectively
[V, T] = sAlign(s1, s2, S, fScore);

% Then, retrieve the alignment steps and the value of the best alignment
steps = retrieveAlignmentSteps(n, m, T);
tscore = V(n+1, m+1);

% Make alignment diagram
ind = cumsum(steps, 2); % convert steps to indices
ind(~ind) = 1; % avoid zero indexing at the begining
aa = [s1(ind(1,:)); s2(ind(2,:))];
if ischar(aa)
    aa(~steps) = '-';
elseif isstring(aa)
    aa(~steps) = "";
elseif isnumeric(aa)
    aa(~steps) = NaN;
end

end

function v = scoreFunc(s1, s2)
% Default scoring function

if isempty(s1) || isempty(s2)
    % Gap
    v = -1;
elseif s1 == s2
    % Match
    v = 1;
%     if s1 == ' '
%         v = 1.5; % higher score for matching spaces
%     end
else
    % Mismatch
    v = -1;
end

end

function [V, T] = sAlign(seq1, seq2, S, scoreFunc)
% this is the main function

n = length(seq1);
m = length(seq2);
pGap = -0.5;

% V stores the values of the best alignments
%   T(i,j) is the value of the best alignment of seq1(1:i-1) and seq2(1:j-1) 
%   under scoring function scoreFunc
V = ones(n+1, m+1) * -inf;

% T stores the best steps; steps are two-dimensional 
T = ones(n+1, m+1, 2) * -inf;

% G stores whether the current best alignment ended with a gap
G = zeros(size(V));

% Fill the table
for i = 1:n+1
    for j = 1:m+1
        % No need to compute the first cell
        if i==1 && j==1
            V(i,j) = 0;
            G(i,j) = 1;
            continue
        end
        
        % Find the best step
        bestVal = -inf;
        bestStep = [-inf -inf];
        for k = 1 : size(S,1)
            % Choose a step
            step = S(k,:);
            a = step(1);
            b = step(2);
            if i-a < 1 || j-b < 1
                % skip out-of-bound step
                continue
            end
            
            % Get the score where we step from
            val = V(i-a,j-b);
            
            % Add penalty if previous step is not gap
            isPrevGap = G(i-a,j-b);
            isThisGap = any(~step); % [0 1] or [1 0]
            if ~isPrevGap && isThisGap
                val = val + pGap;
            end
            
            if val > -inf
                sub1 = seq1(i-a:i-1);
                sub2 = seq2(j-b:j-1);
                scor = scoreFunc(sub1, sub2);
                if val+scor > bestVal
                    bestVal = val+scor;
                    bestStep = [a b];
                end
            end
        end
        
        % Save variables
        V(i,j) = bestVal;
        T(i,j,:) = bestStep;
        if any(~bestStep)
            G(i,j) = 1; % indicate gap
        end
    end
end

end


function ret = retrieveAlignmentSteps(n,m,T)
% retrieve the optimal alignment steps from T

a = T(n+1, m+1, 1);
b = T(n+1, m+1, 2);

if a > -inf && b > -inf
    ret = [retrieveAlignmentSteps(n-a,m-b,T) [a;b]];
else
    ret = [];
end

end

% Original documentation
% 
% Determine the best alignment, with allowed 'steps' <S>, between <sequence1> and <sequence2>
% under similarity function <similarityFunHandle>
% This can equivalently be used for computing a generalized edit distance between two sequences. 
%
% INPUT ARGUMENTS
% <sequence1>: a sequence, e.g., [1,2,2,4] or ['a','b','c','d']
% <sequence2>: another sequence
% <steps>: the allowed steps. If S=[1 0; 0 1; 1 1], the algorithm is identical to the classical 
%                                         Needleman-Wunsch sequence alignment procedure
%            steps must have dimension [K,2] for some K>0.
% <similarityFunHandle>: a function handle (e.g., @myFun) that returns a real number 
%                                                for each of the K steps in S.
%
% RETURN VALUES
% The return value <alignment> is a 2xR matrix with the optimal steps in it.
% For example, the return value
% 1 2 2 3
% 2 1 0 1
% represents the alignment
%    <sequence1> =      x1           x2 x3     x4 x5    x6 x7 x8
%    <sequence2> =    y1 y2           y3                        y4
%
% the return value <value> is the value of the best alignment
%
% Running time of the algorithm is O(nm|S|)
% where n = length(<sequence1>)  and m = length(<sequence2>)
% 
%
% Matlab implementation:
% Steffen Eger, steffen.eger@yahoo.com, 12/3/2013
%
% This implementation is based on 
% Steffen Eger, Sequence alignment with arbitrary steps and further generalizations, with applications to alignments in linguistics. Information Sciences (2013), 237: 287--304.
% see also:
% B. John Oommen, String Alignment With Substitution, Insertion, Deletion, Squashing, and Expansion Operations. Information Sciences (1995), 83: 89--107.