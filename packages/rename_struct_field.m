function s = rename_struct_field(a, oldField, newField)
% function s = rename_struct_field(a, oldField, newField)
% rename a field name inside a struct
[a.(newField)] = a.(oldField);
a = rmfield(a,oldField);
