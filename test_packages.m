clear; close all; clc;

package_list = {'+geo', '+ode', '+opt'};
for p_i = 1:length(package_list)
    package_name = package_list{p_i};
    file_list = dir(package_list{p_i});
    for i = 1:length(file_list)
        if file_list(i).isdir
            continue;
        end

        filename = file_list(i).name;
        idx = regexp(filename, '^test_');
        if ~isempty(idx)
            fprintf('---- ');
            cmd = sprintf('%s.%s()', package_name(2:end), filename(1:end - 2));
            eval(cmd);
        end
    end
end
