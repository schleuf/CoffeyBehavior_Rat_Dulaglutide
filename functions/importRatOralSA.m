function [varTable, eventCode, eventTime] = importRatOralSA(filename)
%IMPORTFILE Import MED-PC data from a text file

%% Input:
%filename = Full file path of Med-PC file generated from Golden_Liq_SelfAdmin.mpc

%% Output
% varTable = table containing variables calulated directly in MED-PC
% including Subject,Session,Date,FileName,EarnedInfusions,TotalInfusions,HeadEntries,...
% Latency,ActiveLever,InactiveLever,Weight
%
% MEP-Array Locations
% \   A(14) = Animal weight (g)
% \   A(15) = Training Session
% \   C(0)  = Infusions
% \   C(4)  = Count of head entries
% \   C(9)  = Current average latency
% \   C(10) = Total active Lev presses
% \   C(11) = Total inactive Lev presses
% \   C(12) = Number of non-contingent reward deliveries
%
% eventCode = raw event code array
%
% eventTime = raw event time array
%
% \  Event stamps ( array E )
% \   1   = Right Lev press
% \   2   = Leftt Lev press
% \   3   = Rewarded rt Lev press
% \   4   = Rewarded lt Lev press
% \   5   = Reinforcement delivery
% \   6   = Head entry
% \   7   = Right Lev light on
% \   8   = Left Lev light on
% \   9   = Right Lev light off
% \   10  = Left Lev ight off
% \   11  = Houselight on
% \   12  = Houselight off
% \   13  = Tone on
% \   14  = Tone off
% \   15  = Magazine tray light on
% \   16  = Magazine tray light off
% \   17  = Infusion turns on
% \   18  = Infusion turns off
% \   19  = Timeout ends and/or drug available at START
% \   20  = Rt Lev press during timeout
% \   21  = Lt Lev press during timeout
% \   22  = Active Lev press
% \   23  = Inactive Lev press
% \   24  = Timeout starts
% \   25  = Time-limit to respond reached, non-contingent reward delivery
% \   26  = Right lever extends
% \   27  = Left lever extends
% \   28  = Right lever retracts
% \   29  = Left lever retracts
% \   30  = Opto Pulse Onset
% \   100 = Session termination

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = [":", " ", "  ", "   ", "    ", "     ", "      ", "       ", "        ", ";"];

%% Import the data
raw = readtable(filename, opts);

%% Find all variables automatically calculated by MED-PC
Date = datetime(raw.Var4(find(raw.Var1=="Start",1,'first')),'InputFormat','MM/dd/yy');
Subject = categorical(raw.Var3(find(raw.Var1=="Subject",1,'first')));
Weight = str2num(raw.Var8{find(raw.Var1=="A")+3}); % Animal Weight 
Session = str2num(raw.Var4{find(raw.Var1=="A")+4}); % Session 

TotalInfusions = str2num(raw.Var4{find(raw.Var1=="C")+1});
EarnedInfusions = TotalInfusions-str2num(raw.Var6{find(raw.Var1=="C")+3});
HeadEntries = str2num(raw.Var8{find(raw.Var1=="C")+1});
Latency = str2num(raw.Var8{find(raw.Var1=="C")+2});
ActiveLever = str2num(raw.Var4{find(raw.Var1=="C")+3});
InactiveLever = str2num(raw.Var5{find(raw.Var1=="C")+3});
FileName = {filename};

varTable = table(Subject,Session,Date,FileName,EarnedInfusions,TotalInfusions,HeadEntries,Latency,ActiveLever,InactiveLever,Weight);

%% Find Raw Event Array and Raw Timestamp Array
eIDX=find(raw.Var1=="E");
tIDX=find(raw.Var1=="T");
e=str2double(reshape(raw{eIDX+1:tIDX-1,4:8}',[],1));
t=str2double(reshape(raw{tIDX+1:end,4:8}',[],1));
eventCode=e(~isnan(e));
eventTime=t(~isnan(t));
end