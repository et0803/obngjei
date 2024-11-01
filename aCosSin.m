function [alpha] = aCosSin(cosValue,sinValue)
% alpha的返回值范围是[-pi,pi)
if sinValue>0
    if cosValue >0
        alpha = asin(sinValue);
    else
        alpha = pi - asin(sinValue);
    end
else
    if cosValue >0
        alpha = asin(sinValue);
    else
        alpha = -pi-asin(sinValue);
    end
end
end
