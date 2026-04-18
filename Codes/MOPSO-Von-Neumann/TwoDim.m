function neighbors=TwoDim(rows, cols, nPop,i)


linear_index =i;

[row, col] = ind2sub([rows, cols], linear_index);

neighbors = [];

if row > 1
    neighbors = [neighbors, sub2ind([rows, cols], row-1, col)];
end

if row < rows
    neighbors = [neighbors, sub2ind([rows, cols], row+1, col)];
end

if col > 1
    neighbors = [neighbors, sub2ind([rows, cols], row, col-1)];
end

if col < cols
    neighbors = [neighbors, sub2ind([rows, cols], row, col+1)];
end


