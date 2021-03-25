function vect = nonNaN2vect(X)

x = permute(X,[3 2 1]);
vect = x(~isnan(x(:)));

