function label = correl_Label(score)
    if score == -1
        label = 'perfect-';
    elseif score < -0.9
        label = 'very strong-';
    elseif score < -0.6
        label = 'strong-';
    elseif score < -0.3
        label = 'moderate-';
    elseif score < 0-eps
        label = 'weak-';
    elseif (score > 0-eps) && (score < 0+eps) 
        label = 'no relationship';
    elseif score == 1
        label = 'perfect+';
    elseif score > 0.9
        label = 'very strong+'  ;  
    elseif score > 0.6
        label = 'strong+';
    elseif score > 0.3
        label = 'moderate+'  ;
    elseif score > 0+eps
        label = 'weak+' ;
    else
        error('Could not find label, check value')
    end
end