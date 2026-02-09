function randvals = getRandomValsFromArray(arrayforrand, Nrand)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
irand    = randperm(numel(arrayforrand), min(numel(arrayforrand), Nrand));
randvals = arrayforrand(irand);
end