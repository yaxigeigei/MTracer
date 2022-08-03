function [ r1, r2, m1, m2 ] = XIndexing( v1, v2 )
%CROSSINDEXING Summary of this function goes here
%   Detailed explanation goes here

% Length of vectors
l1 = length(v1);
l2 = length(v2);

% Length of resulting vector
lr = l1 + l2 - 1;

% Cross-indexing masks for vector 1 and 2
m1 = zeros(lr, l1);
for i = 1 : l1
    m1(i : i+l2-1, i) = 1;
end

m2 = zeros(lr, l2);
for i = 1 : l2
    m2(i : i+l1-1, l2-i+1) = 1;
end

% Cross-indexing range matrices for vector 1 and 2
r1 = zeros(lr, 2);
r2 = zeros(lr, 2);
for i = 1 : lr
    ind = find(m1(i,:));
    r1(i, :) = [ min(ind), max(ind) ];
    ind = find(m2(i,:));
    r2(i, :) = [ min(ind), max(ind) ];
end


end

