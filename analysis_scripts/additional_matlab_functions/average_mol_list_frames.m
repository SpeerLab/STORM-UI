function M2 = average_mol_list_frames(M)

% syntax: M2 = average_mol_list_frames(M)
%
% This function averages the fields in an insight3 style molecule list so
% that there is only one element per molecule. The first element of the
% input fields is omitted (all fields must have at least two elements).


error_text = 'unsuitable input struct M to average_mol_list_frames subroutine';

try
    field_names = fields(M);
    for k = 1:length(M)
        for m = 1:length(field_names)
            M2(k).(field_names{m}) = double(mean(M(k).(field_names{m})(2:end)));
        end
    end
catch
    disp(error_text1)
end
