function generate_rx_power(rec_power)
fileID = fopen('rx_power.txt','w');
num_scenario = size(rec_power,1);
num_combo = size(rec_power,2);
        
for i = 1:num_scenario
    for j = 1:num_combo
        fprintf(fileID,'%.6e,',rec_power(i,j));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);