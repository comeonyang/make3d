function theta = angle(q)
% ANGLE or phase of quaternion.
% (Quaternion overloading of standard Matlab function.)
%
% If q = A + mu .* B where A and B are real/complex, and mu is a unit pure
% quaternion, then angle(q) = B. For pure quaternions, angle(q) = pi/2.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Modified : 20 September 2005 to handle complexified quaternions.

error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

if ispure(q)
    
    % q is a pure quaternion, with no scalar part, and by definition, theqq
    % angle is a right-angle, even for complexified quaternions. Return a
    % matrix of the same size as q, with pi/2 at every position.
    % Rationale : we can construct an exact result in this case with
    % minimal calculation, whereas if we use the algorithm below for full
    % quaternions, the result may be subject to small errors.
    
    theta = ones(size(q)) .* (pi ./ 2);
    
else
    
    % q is a full quaternion. We adopt two methods here, depending on
    % whether q is a real quaternion or a complexified quaternion. In
    % both cases we need the scalar part and the modulus of the vector
    % part. We extract these from unit(q) because the argument or angle
    % of unit(q) is the same as the argument of q, and having unit
    % modulus simplifies the formula in the complex case.

    u =  unit(q);
    x =     s(u);
    y = abs(v(u));

    % Now establish whether x and y are real or complex (strictly whether
    % any imaginary part of x or y is non-zero).
 
    if isreal(x) & isreal(y)

        % All elements of x and y are real. Therefore we can use the
        % simple method of constructing a complex value from x and y
        % and using the Matlab angle function to compute the angle.
        % Because we took the absolute value (modulus) of the vector
        % part when we constructed y, y is never negative, and the
        % angle returned is always between 0 and pi inclusive. This is
        % normal for quaternions, since the vector part is a directionn
        % in 3-space. The opposite direction is represented by negating
        % mu (the axis of the quaternion), not by adding pi to the angle.

        theta = angle(complex(x, y));
    else

        % At least one element of x and/or y is complex, so we must use
        % an algorithm suitable for the complex case. See the appendix
        % below for details of the derivation of the formula used here.
       
        theta = - i .* log(x + i.*y);
    end
end

% Appendix: calculation of arctangent(y, x) when y and x are complex.
%
% Algorithms for computing the arctangent (or arcsin/arccos) of a complex
% value are given in some books on complex analysis. For example:
%
% John H. Mathews,
% 'Basic Complex Variables for Mathematics and Engineering',
% Allyn and Bacon, Boston, 1982. ISBN 0-205-07170-8.
%
% The algorithm used here can be derived using the method given in Mathews'
% book, as follows. We are given x and y, both complex, and these are the
% sine and cosine of theta (also complex), respectively.
% Therefore, tan(theta) = y/x. The sine and cosine of a complex angle can
% be written in terms of complex exponentials as (using mathematical
% notation with implicit multiplication, rather than Matlab notation):
%
% sin(theta) = (1/2i)[exp(i theta) - exp(-i theta)]
% cos(theta) = (1/2) [exp(i theta) + exp(-i theta)]
%
% Therefore, we can write:
%
% y/x = sin(theta)/cos(theta)
%     = [exp(i theta) - exp(-i theta)]/i[exp(i theta) + exp(-i theta)]
%
% Multiplying numerator and denominator on the right by exp(i theta)
% this simplifies to:
%
% y/x = [exp(2i theta) - 1]/i[exp(2i theta) + 1]
%
% and cross multiplying and gathering all terms on the left, and
% multiplying through by -1, we get:
%
% => exp(2i theta)[x - iy] - [x + iy] = 0
%
% => exp(2i theta) = [x + iy]/[x - iy]
%
% Rationalizing the right-hand side we obtain:
%
% exp(2i theta) = [x+iy]^2 / [x^2 + y^2]
%
% Taking the square root of both sides:
%
% exp(i theta) = [x + iy]/sqrt[x^2 + y^2]
%
% Finally, taking the natural logarithm of both sides and multiplying
% throughout by i, we get:
%
% theta = -i ln[(x + iy)/sqrt(x^2 + y^2)]
%
% If we normalise the quaternion before computing x and y, we can ensure
% that sqrt(x^2 + y^2) = 1, and therefore we can simplify the formula to
% that used in the code above.
