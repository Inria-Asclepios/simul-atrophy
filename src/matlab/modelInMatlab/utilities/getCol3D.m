% getCol

function col = getCol3D(i,j,k,m,n)

col = i + (j-1)*m + (k-1)*m*n;

end