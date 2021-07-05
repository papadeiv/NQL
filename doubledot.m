function inner_product = doubledot(vector1,vector2)
    
    vector = reshape(vector1, [], 1).*reshape(vector2, [], 1);
    sorted_vector = sort(vector);
    inner_product = sum(sorted_vector);
end