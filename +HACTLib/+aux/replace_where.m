function arr_out = replace_where(mask, A, B)
	A(~isfinite(A)) = 0;
	B(~isfinite(B)) = 0;

	arr_out = mask .* A + (~mask) .* B;
end