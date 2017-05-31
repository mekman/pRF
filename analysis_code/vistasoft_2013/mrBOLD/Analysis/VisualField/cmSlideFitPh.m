function [distanceShift] = cmSlideFitPh(cortMag)
nMeridia = length(cortMag.corticalDist);
% These are the values that we will return
% These are the shifts we will test.  We look at the total range of the
if(cortMag.templateRoiNum==0 | isempty(cortMag.corticalDist{cortMag.templateRoiNum}))
uDistTemplate = cortMag.corticalDist{cortMag.templateRoiNum};

%mnPhTemplate = unwrapPhases(complexPh2PositiveRad(mnPhTemplate));
        % Set the shift range to be half of the distance range of the data.
        uDist = cortMag.corticalDist{ii};
        %mnPh = unwrapPhases(complexPh2PositiveRad(mnPh));
        % selectGraphWin; plot(uDist,angle(mnPh),'ro')
        
        % Test exhaustively across the shift range
        % 
        err = zeros(size(dShift));
        for jj = 1:length(dShift)
            err(jj) = CMFfitOneMeridianFun(dShift(jj), uDistTemplate, mnPhTemplate, uDist, mnPh);
        end
        
        % selectGraphWin; plot(dShift,err)
        [val idx] = min(err);
        if(min(err)>100)
            % really bad estimate- possibly no overlapping values. Default to no shift.
            distanceShift(ii) = 0;
            disp(['ROI ',num2str(ii),': Distance shift failed- maybe delete this ROI?']);
        else
            distanceShift(ii) = dShift(idx);
            fprintf('ROI %.0f: Distance shift was %.0f mm (range %.0f mm).\n',ii,dShift(idx),r);
        end
    end
end

return;