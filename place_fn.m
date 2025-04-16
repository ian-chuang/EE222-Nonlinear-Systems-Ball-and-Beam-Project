function K = place_fn(A, B, p)
    K = zeros(4, 2);
    K(:) = place(A, B, p);
end