function A = random_transition_matrix(m, n, nonzero_density)

	A = sprand(m, n, nonzero_density);
	A = triu(A) + tril(A');
	A = 100 * A .^ 3;
	A_diag = spdiags(A, 0);
	A = A - A_diag;
	A_offdiag_sums = sum(A-spdiags(A_diag, size(A,1), size(A,2)), 2);
	A = A - spdiags(A_offdiag_sums, 0, size(A,1), size(A,2));
end