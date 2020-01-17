function [c,ceq] = mycon(z)

c   = 0;
ceq = z(2) - z(1).^4 - 2*z(1).^3 + 1.2*z(1).^2 + 2*z(1);

end

