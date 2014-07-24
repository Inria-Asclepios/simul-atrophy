% isCSF, to check whether the position i,j is in CSF or not


function true_false = isCSF(i,j,m,n,csf_w)

if (i < csf_w || j < csf_w || i > (m-csf_w) || j > n-csf_w)
    true_false = true;
else
    true_false = false;
end

end



