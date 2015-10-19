function result = is_in_str(pattern_re, str)
% function result = is_in_str(pattern_re, str)
% returns true is string pattern is in str
% string pattern can be regex expression
result = ~isempty(regexp(str, pattern_re));
end
