
a1 = 1;
b1 = 0.5;

x = -a1/2:0.005:a1/2;
y = -b1/2:0.005:b1/2;

[x, y] = meshgrid(x, y);

a2 = 2;
b2 = 1;

m = 1;
n = 2;

p = 2;
r = 3;

I1_ = (m .* pi./a1) * cos(m.*pi./a1 .* (x + a1./2)) .* cos(p.*pi./a2 .* (x + a2./2))...
    .* (-r.*pi./b2) .* sin(y + b2./2) .* sin(n.*pi/b1 .* (y + b1./2));

I1 = sum(sum(I1_));

I2_ = (n .* pi./b1) * cos(n.*pi./b1 .* (y + b1./2)) .* cos(r.*pi./b2 .* (y + b2./2))...
    .* (-p.*pi./a2) .* sin(x + a2./2) .* sin(m.*pi/a1 .* (x + a1./2));

I2 = sum(sum(I2_));

I = I1 - I2