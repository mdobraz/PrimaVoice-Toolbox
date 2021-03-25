function [ci1,ci2,med,L,U] = CI_values(xpd,csypd,confidence_level)

p1 = (1 - confidence_level) / 2;
p2 = 1 - p1;

[~,x1] = min(abs(csypd - p1));
[~,xm] = min(abs(csypd - 0.5));
[~,x2] = min(abs(csypd - p2));

ci1 = xpd(x1);
med = xpd(xm);
ci2 = xpd(x2);

L = med - ci1; % for usage with funtion errorbar
U = ci2 - med; % for usage with funtion errorbar