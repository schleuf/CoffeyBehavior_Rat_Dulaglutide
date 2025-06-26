function [eventCode,eventTime] = eventFilterHE(eventCode,eventTime)
    % Filters Oral SA events to allow for photometry analysis
    % Removes duplicate events (rapid levers, head entires, etc.)
    % Determines if head entries are rewarded or unrewarded
    
    % INPUT EVENT CODES
    % Event Code 13 = Infusion, Tone, Light, Etc (No need for filtering)
    % Event Code 22 = Active Lever Press
    % Rewarded Press 3 = Active Lever Press Followed by Tone etc.
    % Event Code 23 = Inactive Lever Press (Filter out presses within 10s)
    % Event Code 6 = Head Entry (Filter withing 10s & Rewarded vs Unrewarded)
    
    % OUTPUT EVENT CODES
    % Event Code 13 = Infusion, Rewarded Press, Tone, Light, Etc (No need for filtering)
    % Rewarded Press 3 = Active Lever Press Followed by Tone etc.
    % Event Code 97 = Unrewarded Active Lever Press (Filter out presses within 10s)
    % Event Code 23 = Inactive Lever Press (Filter out presses within 10s)
    % Event Code 96 = Head Entry Filtered Within 2s %SSnote: this is set to 5s
    % Event Code 98 = Rewarded Head Entry (Filter withing 10s)
    % Event Code 99 = Unrewarded Head Entry (Filter withing 10s)
    
    eventTime = round(eventTime,1);
    
    rewLP=eventTime(eventCode==3); 
    cue=eventTime(eventCode==13); 
    actLP=eventTime(eventCode==22);
    inactLP=eventTime(eventCode==23);
    HE=eventTime(eventCode==6);
    rewHE=[];
    unrewHE=[];
    unrewLP=[];
    inLP=[];
    filtHE=[];
    
    % Get rewHE (the first head entry time after each light-on cue)
    c=0;
    curHE=0;
    for i=1:height(cue)
        if ~isempty(HE(find(HE>cue(i),1,'first'))) & HE(find(HE>cue(i),1,'first'))~=curHE
            c=c+1;    
            rewHE(c,1)=HE(find(HE>cue(i),1,'first'));
            curHE=HE(find(HE>cue(i),1,'first'));
        end
    end
    
    % remove rewarded head entry events from HE
    HE = setdiff(HE,rewHE); 
    
    % remove head entry events that are less than 5 seconds before or after a rewarded head entry
    for i=1:height(HE)
        idx(i,1)=~any(abs(rewHE-HE(i)) <= 5); % SSnote: this doesn't seem like it should be absolute value
    end
    HE=HE(idx);  
    
    c=0;
    for i=2:height(HE)
        if ~sum(HE(i)==rewHE) % if HE(i) not in rewHE (SSnote: how could it be...? HE is set to all non rewHE values of previous HE var)
            if HE(i)-HE(i-1)>5
                c=c+1;
                unrewHE(c,1)=HE(i);
            end
        end
    end

    
    % if active lever press is unrewarded and more than 10 seconds after the 
    % previous active lever press, put it in the unrewarded active lever press bin
    % SSnote: why though? shouldn't all active presses ten seconds apart be
    % rewarded? doesn't this drop a lot of unrewarded lever presses?
    c=0;
    for i=2:height(actLP)
        if ~sum(actLP(i)==rewLP) 
            if actLP(i)-actLP(i-1)>10 
                c=c+1;
                unrewLP(c,1)=actLP(i);
            end
        end
    end
    
    % discount inactive lever presses within 10 seconds of eachother
    if ~isempty(inactLP)
        c=1;
        inLP(c,1)=inactLP(1);
        for i=2:height(inactLP)
            if inactLP(i)-inactLP(i-1)>10
                c=c+1;
                inLP(c,1)=inactLP(i);
            end
        end
    end
    
    % filtHE = unrewarded head entries that are at least 5 seconds apart 
    % from each other
    % SSnote: why not just do this first before extracting rewarded head entries? 
    % don't we want rewarded head entries included in filtered head entries??
    % isn't this the same as unrewHE
    if ~isempty(HE)
        c=1;
        C(c,1)=HE(1);
        for i=2:height(HE)
            if HE(i)-HE(i-1)>5
                c=c+1;
                filtHE(c,1)=HE(i);
            end
        end
    end
    
    eventCode=[repmat(3,[height(rewLP),1]);... 
        repmat(13,[height(cue),1]);...
        repmat(97,[height(unrewLP),1]);...
        repmat(23,[height(inLP),1]);...
        repmat(96,[height(filtHE),1]);...
        repmat(98,[height(rewHE),1]);...
        repmat(99,[height(unrewHE),1])];
    
    eventTime=[rewLP;... % rewarded lever presses with head entry before next rewarded lever press
        cue;... % all rewarded lever presses
        unrewLP;... % unrewarded lever presses > 5s after a previous active press
        inLP;... % inactive lever presses > 10s apart from previous inactive press
        filtHE;... % unrewarded head entries that are at least 5s apart from each other
        rewHE;... % head entries immediately following a rewarded lever press
        unrewHE]; % redundant with fileHE
end