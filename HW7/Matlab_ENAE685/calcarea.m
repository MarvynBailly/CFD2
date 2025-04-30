% calcarea.m
%
% INPUTS: x (location of points along the nozzle, array)
% 
% OUTPUTS: area (area at each of the locations, array)
%
function area = calcarea(x)
   area = 1.398+0.347*tanh(0.8*(x-4.0));
end
