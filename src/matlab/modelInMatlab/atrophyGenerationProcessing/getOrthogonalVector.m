
function [out_v] = getOrthogonalVector(v)

if(norm(v) < eps)
    out_v = v;
    return
elseif(v(1)>eps || v(2)>eps)
    out_v(1) = v(2);
    out_v(2) = -v(1);

else
    out_v(1) = 1;
    out_v(2) = 0;
end

out_v(3) = 0;
out_v = out_v./norm(out_v);
out_v(out_v==Inf) = 0;