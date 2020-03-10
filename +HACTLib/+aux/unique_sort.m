function [values_out, cdf_out, iu] = unique_sort(values_in, pmf_in, iunique)
	sorted_mat = sortrows([values_in(:), pmf_in(:)]);
	tmp_cdf = cumsum(sorted_mat(:,2));

	if iunique == 1
		[values_out, iu] = unique(sorted_mat, 'last');
		cdf_out = tmp_cdf(iu);
	else
		[cdf_out, iu] = unique(tmp_cdf, 'first');
		values_out = sorted_mat(iu,1);
	end
end