% Better making recursive, for now just a loop used instead of recursion!
% computes exponential map of a field v: exp(v) by using v/(2^n) as the
% smallest displacement that allows us to approximate adding of
% displacements succesively to follow a material point.

function vd = expRecursive(v,n)

vd = v/(2^n);
for indx = 1:n
    vd = vecCompose(vd,vd);
end
end

