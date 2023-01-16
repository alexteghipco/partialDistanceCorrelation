function [tmpl] = getord(tmpo)

counter = 1;
clear tmpl
for j = 1:length(tmpo)
    if ~isempty(tmpo)
        if rem(counter,2) == 1
            tmpl(j) = tmpo(end);
            tmpo(end) = [];
        else
            tmpl(j) = tmpo(1);
            tmpo(1) = [];
        end
        counter = counter+1;
    end
end