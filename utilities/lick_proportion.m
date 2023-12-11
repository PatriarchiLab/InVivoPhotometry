function meanLickProportion = lick_proportion(lick_vectors_cells_sound,lick_vectors_cells_No_sound)

for ii=1:size(lick_vectors_cells_sound,1)
   
    Lick_Sound = nonzeros(lick_vectors_cells_sound(ii,:));
    LickPercentSound(ii)=100*(length(Lick_Sound)/size(lick_vectors_cells_sound,2));
    
    Lick_NoSound = nonzeros(lick_vectors_cells_No_sound{1,ii}); 
    if length(lick_vectors_cells_No_sound{1,ii})==0
        LickPercentNoSound(ii)=0;
    else
    LickPercentNoSound(ii)=100*(length(Lick_NoSound)/length(lick_vectors_cells_No_sound{1,ii})); 
    end

end

meanLickProportion=(mean(LickPercentSound))/(mean(LickPercentNoSound));

end