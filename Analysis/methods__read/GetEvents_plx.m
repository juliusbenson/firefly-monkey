function [t_events] = GetEvents_plx(fname)
% get begin, reward, and end times from plx file

[~, ts, sv] = plx_event_ts(fname, 257); count=0;
t_events.start = ts(sv==1);
ts(sv==1) = []; sv(sv==1) = [];
ts(sv==8) = []; sv(sv==8) = [];
nevents = find(sv==3,1,'last');
for i=1:nevents
    if sv(i)==3
        count=count+1;
        t_end(count)=ts(i);
        if sv(i-1)==4
            t_rew(count)=ts(i-1);
            t_beg(count)=ts(i-2);
        elseif sv(i-1)==2
            t_rew(count)=nan;
            t_beg(count)=ts(i-1);
        elseif sv(i-1)==3
            t_rew(count)=nan;
            t_beg(count)=ts(i-1);
        else
            %error('strobed unknown entity');
        end
    end
end

t_events.t_beg = t_beg;
t_events.t_rew = t_rew;
t_events.t_end = t_end;

% convert eventtimes from samples to seconds 
%t_events.start = t_events.start/fs;
%t_events.t_beg = t_events.t_beg/fs; t_events.t_end = t_events.t_end/fs; t_events.t_rew = t_events.t_rew/fs;
t_events.start = t_events.start;
t_events.t_beg = t_events.t_beg; t_events.t_end = t_events.t_end; t_events.t_rew = t_events.t_rew;