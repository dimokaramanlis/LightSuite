function A = copyStructBtoA(A, B)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for fn = fieldnames(B)'
    % if isfield(A, fn)
    %     warning('found common stuct names');
    % end
    A.(fn{1}) = B.(fn{1});
end

end