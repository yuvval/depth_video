
x = dir('*.m'); 
[info,files]=mlint({x.name},'-fullpath');
for i=1:length(info)
  if ~isempty(info{i})
    fprintf('\n===============%s:\n',files{i});        
    ml_temp_data = info{i};
    msgs = {ml_temp_data.message};
    lines = {ml_temp_data.line};
    cols = {ml_temp_data.column};
    n_msgs = length({ml_temp_data.message});

    for j=1:n_msgs
      fprintf('line=%3d ',lines{j});
      % fprintf('col=%d',cols{j});
      fprintf('%s\n',msgs{j});
    end
  end
end
clear ml_temp_data