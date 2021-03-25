function cursor_pos = get_cursor_pos(sl_slice)

sl_fields = fieldnames(sl_slice);
cursor_pos = nan(numel(sl_fields),1);
for i = 1:numel(sl_fields)
    cursor_pos(i) = sl_slice.(sl_fields{i}).Value;
end

cursor_pos = round(cursor_pos);