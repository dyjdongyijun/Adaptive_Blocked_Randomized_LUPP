function d = square_difference(x, y)
    w = x + y;
    z = x - y;

    z(z<0) = 0;
    d = sqrt(w .* z);
end