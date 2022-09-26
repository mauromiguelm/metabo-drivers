columns = 1:24

Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
rows_numbers = [1:16];
rows = Alphabet(rows_numbers)

clear Alphabet numbers

layout = matAssayLayouts(:,:,1);

store_matrix = cell(384, length(cellHeader)+2)

count = 0

for row = rows_numbers 
        
        for col = columns
            
            count = count +1;
            
            %row = 1;
            %col = 1;
            
            store_matrix{count,1} = rows(row);
            store_matrix{count,2} = col;
            %store_matrix{count,3:size(store_matrix,2)} =  (cellConditions{layout(row,col),:})
            
            if ~isnan(layout(row,col))
             
                store_matrix(count,3:size(store_matrix,2)) =  cellConditions(layout(row,col),:);
                
            end
            
            
                        
        end
end


xlswrite('randomized_layout_2MSP_batch2.xls',store_matrix)
