function output_matrix(f,n,x,y,file)
%output_matrix(f,n,x,y,'field-T-000000.bin'));
% Allows to write as in output_fields() in BASILISK, where
% the file has been stored as a single precision binary file.

fileID = fopen(file,'w');

fwrite(fileID,n,'single');

for j = 1:n
	fwrite(fileID,y(j),'single');
end

for i = 1:n
	fwrite(fileID,x(i),'single');
	for j = 1:n
		fwrite(fileID,f(j,i),'single');
	end
end
