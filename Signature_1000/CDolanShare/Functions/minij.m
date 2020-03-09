function [i,j]=minij(A)
    [Aj,ia]=min(A);
    [~,j]=min(Aj);
    i=ia(j);
    
    return;
end
