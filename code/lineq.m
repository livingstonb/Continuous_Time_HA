function y = lineq(x1, x2, y1, y2, x)
    m = (y2 - y1) / (x2 - x1);
    b = y1 - m * x1;
    y = b + m * x;
end