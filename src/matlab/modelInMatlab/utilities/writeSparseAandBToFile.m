% Function to write sparse matrix A and rhs b of the corresponding linear
% system Ax = b in two ascii files.
% The top line of A matrix has: m n nnz
% 
% 
% 
% Bishesh Khanal
% Asclepios, INRIA

function writeAandBtoFile(A,b,AFilename,bFilename)

[i_indx j_indx x] = find(A);
i_indx = int32(i_indx);
j_indx = int32(j_indx);

% fid_rhs = fopen('petsC/b.txt','wb');
fid_rhs = fopen(bFilename,'wb');
fprintf(fid_rhs, '%f\n',b);
fclose(fid_rhs);

% fid = fopen('petsC/A.txt', 'wb');
fid = fopen(AFilename, 'wb');
% Write the sizes and number of non-zeros at the top of the file.
fprintf(fid,'%d %d %d \n', int32(size(A,1)), int32(size(A,2)), int32(size(x,1)));
for ii = 1:size(i_indx)
        fprintf(fid, '%d %d %f \n', i_indx(ii), j_indx(ii), x(ii));
end
fclose(fid);

end
