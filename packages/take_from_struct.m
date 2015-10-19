function [value,out_options]=take_from_struct(options,fieldname,default)
%
% [VAL,OUT_OPTIONS] = TAKE_FROM_STRUCT(OPTIONS, FIELDNAME, DEFAULT)  
% 
% Take values from the options structure, use default if fiekd does
% not exist. Provide meaningful error messages. The function
% also updates the structure when the default is used.
% If default is not given, program aborts if field does not exist.
%
% Examples: 
% 
%  1. get the values of n_restarts, use default of 20 if field isnt set
%     [options,n_restarts] = take_from_struct(options,'n_restarts',20);
% 
%  2. get the values of n_restarts, abort if doesn't exist
%     [options,n_restarts] = take_from_struct(options,'n_restarts');
%
% (C) Gal Chechik, 2004
% Software available for academic use. Other uses require explicit permission.
% 
  
  out_options = options;    

  try
    value = out_options.(fieldname);
  catch
    if(exist('default','var')==0)
      fprintf('\n\nError:\n');
      fprintf('    Field "%s" does not exist in structure\n\n', fieldname);
      error('Trying to read from a field that does not exist');
    end    
    value=default;
    out_options.(fieldname) = value;
  end
  
  return
end