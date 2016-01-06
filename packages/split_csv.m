function splitted_cell = split_csv(csv_string)
% function splitted_cell = split_csv(csv_string)
% split a csv string to a cell array of strings
splitted_cell = regexp(csv_string, ',?\s?', 'split');

