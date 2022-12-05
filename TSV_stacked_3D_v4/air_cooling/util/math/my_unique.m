function [new_array, len]= my_unique(array, thre)
%this function is to return a unique new_array of array based on the
%user-defined threshold
    len_old = length(array);
    len = len_old;
    new_array = sort(array);
    max_v = max(array);
    
    for i=1:len_old-1
        if new_array(i+1) - new_array(i) < thre
            new_array(i) = max_v+1; 
            %make sure this value will be larger than the max(array)
            len = len - 1;
        end
    end
    
    new_array = sort(new_array);
    new_array = new_array(1:len);

end

