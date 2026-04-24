function scenario = yaml_read(filepath)
%YAML_READ Simple YAML parser for scenario files

fid = fopen(filepath, 'r');
if fid == -1
    error('Cannot open file: %s', filepath);
end

scenario = struct();
current_parent = '';
multiline_key = '';
multiline_indent = 0;

while ~feof(fid)
    line = fgetl(fid);

    if isempty(line) || startsWith(strtrim(line), '#')
        continue;
    end

    % Remove inline comments
    comment_idx = find(line == '#', 1, 'first');
    if ~isempty(comment_idx)
        line = line(1:comment_idx-1);
    end

    stripped = strtrim(line);
    % Count leading whitespace only. Using length(line)-length(stripped)
    % includes trailing whitespace (e.g. left over after inline-comment
    % stripping), which would mis-classify a root-level `key: val  # c`
    % line as nested.
    leading = regexp(line, '^\s*', 'match', 'once');
    indent = length(leading);

    % Handle multiline strings (|)
    if contains(stripped, ': |')
        key = strtrim(extractBefore(stripped, ':'));
        if indent == 0
            scenario.(key) = '';
            multiline_key = key;
        else
            scenario.(current_parent).(key) = '';
            multiline_key = [current_parent, '.', key];
        end
        multiline_indent = indent;
        continue;
    end

    % Continue multiline string
    if ~isempty(multiline_key) && indent > multiline_indent
        if contains(multiline_key, '.')
            parts = strsplit(multiline_key, '.');
            scenario.(parts{1}).(parts{2}) = [scenario.(parts{1}).(parts{2}), stripped, ' '];
        else
            scenario.(multiline_key) = [scenario.(multiline_key), stripped, ' '];
        end
        continue;
    else
        multiline_key = '';
    end

    if indent == 0
        current_parent = '';
    end

    if contains(stripped, ':')
        colon_idx = find(stripped == ':', 1, 'first');
        key = strtrim(stripped(1:colon_idx-1));
        rest = '';
        if colon_idx < length(stripped)
            rest = strtrim(stripped(colon_idx+1:end));
        end

        if ~isempty(rest)
            value = parse_yaml_value(rest);
            if indent == 0
                scenario.(key) = value;
            elseif indent > 0 && ~isempty(current_parent)
                if ~isfield(scenario, current_parent) || ~isstruct(scenario.(current_parent))
                    scenario.(current_parent) = struct();
                end
                scenario.(current_parent).(key) = value;
            end
        else
            if indent == 0
                scenario.(key) = struct();
                current_parent = key;
            end
        end
    end
end

fclose(fid);
end

function value = parse_yaml_value(str)
%PARSE_YAML_VALUE Parse a YAML value string into MATLAB type

value = str;

if startsWith(value, '"') && endsWith(value, '"')
    value = value(2:end-1);
    return;
end

if startsWith(value, '[') && endsWith(value, ']')
    inner = value(2:end-1);
    inner = strrep(inner, 'null', 'NaN');
    arr_parts = strsplit(inner, ',');
    arr = cell(1, length(arr_parts));
    for i = 1:length(arr_parts)
        v = strtrim(arr_parts{i});
        if startsWith(v, '"') && endsWith(v, '"')
            arr{i} = v(2:end-1);
        else
            num = str2double(v);
            if ~isnan(num) || strcmp(v, 'NaN')
                arr{i} = num;
            else
                arr{i} = v;
            end
        end
    end
    value = arr;
    return;
end

num = str2double(value);
if ~isnan(num)
    value = num;
end
end
