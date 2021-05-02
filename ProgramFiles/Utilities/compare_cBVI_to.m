function compare_cBVI_to(s,p,fig_num, is_square)
    % sanitize
    if ~exist('fig_num', 'var')
        fig_num = 1;
    end
    if ~exist('is_square', 'var')
        is_square = false;
    end
    % compute third-order effects
    [cBVI_to, cBVI_opt_to] = calc_tlb_thirdorder(s,p,is_square);
    
    % compute distance ratios
    [dist, dist_opt] = deal(zeros(2, length(p.cBVI)));
    for shch = 1:length(p.cBVI)
        % get endpoints, path length
        true = p.G_locus{shch}{end}.G(end,:);
        true_opt = p.G_locus{shch}{end}.G_opt(end,:);
        path_length = p.G_locus{shch}{end}.S(end);
        % save orig, opt
        dist(1,shch) = norm(true - p.cBVI{shch})/path_length;
        dist(2,shch) = norm(true - p.cBVI{shch} - cBVI_to{shch})/path_length;
        dist_opt(1,shch) = norm(true_opt - p.cBVI_opt{shch})/path_length;
        dist_opt(2,shch) = norm(true_opt - p.cBVI_opt{shch} - cBVI_opt_to{shch})/path_length;
    end
    figure(fig_num);
    subplot(2,1,1);
    plot(dist'*100, '-o');
    xlabel('Path');
    ylabel('Error (% of length)');
    legend('cBVI', 'cBVI_{to}');
    subplot(2,1,2);
    plot(dist_opt'*100, '-o');
    xlabel('Path');
    ylabel('Error (% of length)');
    legend('cBVI', 'cBVI_{to}');
end