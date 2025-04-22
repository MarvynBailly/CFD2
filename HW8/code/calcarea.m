% calcarea.m

function area = calcarea(x)
   % area = 0.05*(1.+4.*(x-0.5).^2);
   area = 1.398 + 0.347 * tanh(0.8*(x - 4));
end
