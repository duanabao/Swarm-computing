function ys = raicosfct(bws, ts);

r = 0.5;
fo = bws/2;
df = r*fo;

if ts == 0,
    ys = 1;
elseif abs(df*ts) == 0.25,
    ys = (sin(2*pi*fo*ts)/(2*pi*fo*ts)) * (pi/(16*df*ts));
else,
    ys = (sin(2*pi*fo*ts)/(2*pi*fo*ts)) * (cos(2*pi*df*ts)/(1 - (4*df*ts)^2));
end;
return;