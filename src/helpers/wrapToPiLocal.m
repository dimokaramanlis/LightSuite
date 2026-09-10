function ang = wrapToPiLocal(ang)
%WRAPTOPILOCAL Wrap angles to (-pi, pi] without any toolbox dependency.
%
%   ANG = WRAPTOPILOCAL(ANG) is equivalent to wrapToPi (Mapping Toolbox) and to
%   circ_dist(ANG, 0) (CircStat), but relies only on base MATLAB so the spinal
%   cord pipeline runs without either toolbox installed.

ang = angle(exp(1i * ang));

end
