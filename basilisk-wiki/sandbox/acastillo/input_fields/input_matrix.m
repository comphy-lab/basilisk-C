%% input_matrix
% Allows to read the output generated from output_fields() in BASILISK, where
% the file has been stored as a single precision binary file.
%
% [T,n,x,y] = input_matrix('field-T-000000.bin'));
%
function [f,n,x,y] = input_matrix(file)

fileID = fopen(file,'r');
n = fread(fileID,1,'single');
y = zeros(n,1); x = zeros(n,1); f = zeros(n,n);

for j = 1:n
	y(j) = fread(fileID,1,'single');
end

for i = 1:n
	x(i) = fread(fileID,1,'single');
	for j = 1:n
		f(j,i) = fread(fileID,1,'single');
	end
end
