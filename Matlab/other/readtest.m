% количество элементов в ячейке
block_size = 10;

% describe the format of the data
% for more information, see the textscan reference page
format = '%d %d';

file_id = fopen('test.txt');

while ~feof(file_id)
   segarray = textscan(file_id, format, block_size, 'delimiter', '\n');   
end

fclose(file_id);