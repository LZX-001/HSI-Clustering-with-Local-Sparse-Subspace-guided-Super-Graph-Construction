function adj_sp = createGraph(label,row,col)
%CREATEGRAPH 此处显示有关此函数的摘要
%   此处显示详细说明

label_sp_M = reshape(label, [row, col]);
num_sp=max(label(:));
adj_sp_temp = zeros(num_sp);
for ii_sp = 1:num_sp-1
    [rr_1, cc_1] = find(label_sp_M == ii_sp);
    rc_1 = [rr_1, cc_1];
    for jj_sp = ii_sp+1:num_sp
        [rr_2, cc_2] = find(label_sp_M == jj_sp);
        rc_2 = [rr_2, cc_2];
        dist_temp = pdist2(rc_1, rc_2);
        if min(dist_temp(:)) == 1
            adj_sp_temp(ii_sp, jj_sp) = 1;
        end
    end
end

adj_sp = adj_sp_temp + adj_sp_temp';
end

