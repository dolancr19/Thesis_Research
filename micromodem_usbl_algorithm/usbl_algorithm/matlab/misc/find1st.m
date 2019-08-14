function        ps = find1st(dv,thresh)
%
%       ps=find1st(dv,thresh)
%       Find the first peak in vector dv that exceeds a threshold.
%
%       dv      is the input 'decision' sequence.
%       thresh  is the level above which a detection is made.
%
%       ps      is the index of the first peak in dv, or -1
%		if thresh is not exceeded in the data record
%
%
%       Mark Johnson    Acoustic Telemetry Group, WHOI
%
%       Last Modified:  31 November 1995
%

% $Id: find1st.m 10422 2013-07-09 19:37:07Z andrew $

% find first point that exceeds the threshold
ps = min(find(dv>=thresh)) ;

% find the immediately following peak

if length(ps)==1,
	ps = ps + min(find(diff(dv(ps:size(dv,1)))<0)) - 1 ;
	if isempty(ps)
	  ps = -1;
	end
	

else
	ps = -1 ;
end
