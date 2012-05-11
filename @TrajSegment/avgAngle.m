function res = avgAngle(obj,units)
% 1/2 acos( <2 cos(theta) > ) as a function of time
%   where theta are the differences between adjacent rings
% units = radians or degrees

if (nargin < 2)
   units = 'radians';
else
   units = validatestring(units,{'radians','degrees'});
end

cosDiff = cos(2.0 * obj.angleDiffs);
t1 = mean(cosDiff,1);
if (units == 'radians')
   res = 0.5 * acos(t1);
else
   res = (0.5 * 180/pi) * acos(t1);
end