function d = greatCircleDistance(phi_s, lambda_s, phi_f, lambda_f, r)
% compute the great circle distance given lat and long for two points
% optionally, a fifth parameter (r) can be specified. If this paramter
% isn't specified it's assumed to be the mean radius of the earth. The
% calculation is done using the Vincenty formula.
%
% INPUTS:
% phi_s    = latitude of the standpoint (base) [rad]
% lambda_s = longitude of the standpoint (base) [rad]
% phi_f    = latitude of the forepoint (destination) [rad]
% lambda_f = longitude of the forepoint (destination) [rad]
% r        = radius of the sphere [units determine units of d]
%
% OUTPUT:
% d        = great circle distance from standpoint to forepoint
%
% See http://en.wikipedia.org/wiki/Great-circle_distance

% If no arguments, bail out
if nargin < 4
    fprintf('Usage: greatCircleDistance(phi_s, lambda_s, phi_f, lambda_f, r)\n')
    return
end

% If no radius supplied, assume the mean radius of the earth in km
if nargin < 5
    r = 6371.01; % km
end

% Compute Delta lambda (delta longitude)
Delta_lambda = lambda_f - lambda_s;

% Compute Delta sigma (central angle)
Delta_sigma = atan2(sqrt((cos(phi_f).*sin(Delta_lambda)).^2 + (cos(phi_s).*sin(phi_f) - sin(phi_s).*cos(phi_f).*cos(Delta_lambda)).^2), ...
    sin(phi_s).*sin(phi_f) + cos(phi_s).*cos(phi_f).*cos(Delta_lambda));

d = r*Delta_sigma;

end