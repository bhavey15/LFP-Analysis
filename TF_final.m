%Note: This script needs the folders sessionInfo and MatFiles to be added to current path

clear all
clc

'This code gives spectrogram and spectrum of all the three classes of both linear and spiral motion stimulus for 34 neurons'

addpath('C:\Users\ACER\Documents\LFP Analysis\Data\MatFiles')

addpath('C:\Users\ACER\Documents\LFP Analysis\Data\Basic tools to read Dan Data\sessionInfos')

addpath('C:\Users\ACER\Documents\LFP Analysis\Data\Basic tools to read Dan Data')

load Dan_sessionsInfo.mat
load Dan_ValidLFPSessions.mat
load Dan_MTValidElecs.mat
load Dan_PrefTable.mat
load MST_LS.mat

nsessions = length(sessionsInfo); %total num of sessions

nelecs = 9; %Note that there are only 3 electrodes recorded from (numbers 1, 5, 9)

MT_DELAY = 0; %80   %an assumptive variable addressing the delay it takes for sensory information to reach area MT (you can set it 0 if your analysis period starts from stimulus onset)
STIM_ONSET = 1437;   %interval after trial start taking before the onset of RDPs // changed later
MIN_PRES_TIME = 1200; %680

%Classes to look into (in all classes, pref direction is shown inside RF):
%34: outside RF: pref,      attention: inside   RF
%26: outside RF: pref,      attention: outside  RF
%28: outside RF: Antipref,  attention: outside  RF
lookingClassesSpiral = [26 28 34]; 
lookingClassesLinear1 = [50 48 49];
lookingClassesLinear2 = [48 46 47];
%Input start and end sessions (for both monkeys nic and wal)
% 1-37: nic; 38-86: wal
params.Fs=1000; % sampling frequency
params.fpass=[0 200]; % band of frequencies to be kept
params.tapers=[3 5]; % taper parameters
params.pad=0; % pad factor for fft
params.err=[2 0.05];
params.trialave=1;
movingwin=[0.2 0.01];
srate=1000;
inp_ss = input('Enter start session, please: ');
inp_es = input('Enter stop  session, please: ');
for isession = inp_ss:inp_es
    
    sessionID = MST_LS(isession);
    strcat('loading session # ', num2str(isession), '...')
    load(char(strcat(sessionID,'.mat')));
    'loaded!'
    fields=fieldnames(sessionsInfo);
    for session = 1:length(sessionsInfo)
        if strcmp(sessionsInfo(session).sessionID,sessionID)
            break;
        end
    end
    disp(session);
    sessionInfo = sessionsInfo(isession).sessionInfo; %Information on the timing of each trial
    ntrials = length(sessionInfo);
    

    %% %Scanning through electrodes and trials
    for ielec = 1:3 
        
        %exists in Sonia's list of valid elecs?
        if MTValidElecs(isession,ielec)==0 
            continue;
        end
        
        for iunit=1:2 %Note that the unit' index (in mat datafiles) starts from 2 (1 is unsorted data)
            %exists in Sonia's list of valid (MT) neurons? (if not, skip to the next neuron)
            if PrefTable(isession,ielec+(iunit-1)*8) == -1
                continue;
            end
            
            %%% added by Sonia to take LFPs of all trial in a given cell.
            %%% modifications made by Bhavey to take LFPs of all trials for
            %%% three different classes(26, 28, 34)
             trialNum_26 = 0;
             trialNum_28 = 0;
             trialNum_34 = 0;
             trialNum_50 = 0;
             trialNum_48 = 0;
             trialNum_49 = 0;
             
            for itrial = 1:ntrials %trial' index as in trialInfo table
                itrial_lfp = find(sumPlxMatch.matchMLTrialID==itrial); %find the index of matching lfp trial
                if isempty(itrial_lfp) % does the matching plx trial (LFP/spikes) trial exist?
                    continue;
                end
%                 disp(itrial);
                %correct and wanted trial and class and long-enough trial, respectively?
                iclass = sessionInfo(itrial,2);
                

                if ((sessionInfo(itrial,3) == 22) && (~isempty(find(lookingClassesSpiral==iclass, 1))) && (sessionInfo(itrial,7) > STIM_ONSET+MIN_PRES_TIME))
                    %trialNum = trialNum+1;
                    %define start and end time of your analysis period
                    %extract lfp and spikes within the favorite interval
                    %(tstart:tend); also a place you could extract the
                    %spike waveforms
                    tstart = (-sessionInfo(itrial,4)+STIM_ONSET+MT_DELAY);
                    tend = (-sessionInfo(itrial,4)+STIM_ONSET+MIN_PRES_TIME);
                    trial_data=lfps(1,1+(ielec-1)*4).LFP{1,itrial_lfp};
                    %uncomment line below for baseline normalization
%                     [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    %normalization by max. Comment this for baseline
                    %normalization
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %uncomment below lines for baseline normalization of
                    %spectrogram
%                     for i=1:length(f)
%                         baseline_p=P(i);
% %                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    [fft1,freq,err]=mtspectrumc(trial_data(1:1200),params);
                    x = struct('S',S,'t',t,'f',f,'Serr',Serr);
                    %uncomment lines below for baseline normalization of
                    %spectrum
%                     for i=1:length(freq)
%                         factor=interp1(F,P,freq(i));
%                         fft1(i)=10*log10(fft1(i)/factor);
%                     end
                    if iclass==26     %if itrial belongs to class 26
                        trialNum_26 = trialNum_26+1;
                        cur_lfp_26(trialNum_26) = x;
                        pf_26(trialNum_26,1:410)=fft1';
                    elseif iclass==28    % if itrial belongs to class 28
                        trialNum_28 = trialNum_28+1;
                        cur_lfp_28(trialNum_28) = x;
                        pf_28(trialNum_28,1:410)=fft1';
                    elseif iclass==34    % if itrial belongs to class 34
                        trialNum_34 = trialNum_34+1;
                        cur_lfp_34(trialNum_34) = x;
                        pf_34(trialNum_34,1:410)=fft1';
                    end
                end
                if ((sessionInfo(itrial,3) == 22) && (~isempty(find(lookingClassesLinear1==iclass, 1))) && isession<=19 && (sessionInfo(itrial,7) > STIM_ONSET+MIN_PRES_TIME))
                    tstart = (-sessionInfo(itrial,4)+STIM_ONSET+MT_DELAY);
                    tend = (-sessionInfo(itrial,4)+STIM_ONSET+MIN_PRES_TIME);
                    trial_data=lfps(1,1+(ielec-1)*4).LFP{1,itrial_lfp};
                     %uncomment line below for baseline normalization
%                   [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    %normalization by max. Comment this for baseline
                    %normalization
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %uncomment below lines for baseline normalization of
                    %spectrogram
%                     for i=1:length(f)
%                         baseline_p=P(i);
% %                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    [fft1,freq,err]=mtspectrumc(trial_data(1:1200),params);
                    x = struct('S',S,'t',t,'f',f,'Serr',Serr);
                    %uncomment lines below for baseline normalization of
                    %spectrum
%                     for i=1:length(freq)
%                         factor=interp1(F,P,freq(i));
%                         fft1(i)=10*log10(fft1(i)/factor);
%                     end
                    if iclass==50     %if itrial belongs to class 26
                        trialNum_50 = trialNum_50+1;
                        cur_lfp_50(trialNum_50) = x;
                        pf_50(trialNum_50,1:410)=fft1';
                    elseif iclass==48    % if itrial belongs to class 28
                        trialNum_48 = trialNum_48+1;
                        cur_lfp_48(trialNum_48) = x;
                        pf_48(trialNum_48,1:410)=fft1';
                    elseif iclass==49    % if itrial belongs to class 34
                        trialNum_49 = trialNum_49+1;
                        cur_lfp_49(trialNum_49) = x;
                        pf_49(trialNum_49,1:410)=fft1';
                    end
                end
                if ((sessionInfo(itrial,3) == 22) && (~isempty(find(lookingClassesLinear2==iclass, 1))) && isession>19 && (sessionInfo(itrial,7) > STIM_ONSET+MIN_PRES_TIME))
                    tstart = (-sessionInfo(itrial,4)+STIM_ONSET+MT_DELAY);
                    tend = (-sessionInfo(itrial,4)+STIM_ONSET+MIN_PRES_TIME);
                    trial_data=lfps(1,1+(ielec-1)*4).LFP{1,itrial_lfp};
%                    %uncomment line below for baseline normalization
%                     [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    %normalization by max. Comment this for baseline
                    %normalization
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %uncomment below lines for baseline normalization of
                    %spectrogram
%                     for i=1:length(f)
%                         baseline_p=P(i);
% %                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    [fft1,freq,err]=mtspectrumc(trial_data(1:1200),params);
                    x = struct('S',S,'t',t,'f',f,'Serr',Serr);
                    %uncomment lines below for baseline normalization of
                    %spectrum
%                     for i=1:length(freq)
%                         factor=interp1(F,P,freq(i));
%                         fft1(i)=10*log10(fft1(i)/factor);
%                     end
                    if iclass==48     %if itrial belongs to class 26
                        trialNum_50 = trialNum_50+1;
                        cur_lfp_50(trialNum_50) = x;
                        pf_50(trialNum_50,1:410)=fft1';
                    elseif iclass==46    % if itrial belongs to class 28
                        trialNum_48 = trialNum_48+1;
                        cur_lfp_48(trialNum_48) = x;
                        pf_48(trialNum_48,1:410)=fft1';
                    elseif iclass==47    % if itrial belongs to class 34
                        trialNum_49 = trialNum_49+1;
                        cur_lfp_49(trialNum_49) = x;
                        pf_49(trialNum_49,1:410)=fft1';
                    end
                end
            end
        end
    end
    %spectrum average of all 6 classes: 3 for LMS and 3 for SMS
    if(exist('pf_26','var'))
        pf_26_avg(isession,1:410)=mean(pf_26);
    end
    if(exist('pf_28','var'))
        pf_28_avg(isession,1:410)=mean(pf_28);
    end
    if(exist('pf_34','var'))
        pf_34_avg(isession,1:410)=mean(pf_34);
    end
    if(exist('pf_48','var'))
        pf_48_avg(isession,1:410)=mean(pf_48);
    end
    if(exist('pf_49','var'))
        pf_49_avg(isession,1:410)=mean(pf_49);
    end
    if(exist('pf_50','var'))
        x=size(pf_50);
        if(x(1)==1)
            pf_50_avg(isession,1:410)=pf_50;
        else
            pf_50_avg(isession,1:410)=mean(pf_50);
        end
    end
    %spectrogram average for all 6 classes: 3 LMS and 3 SMS
    if(exist('cur_lfp_28','var'))
    per_lfp_28(isession) = average(cur_lfp_28);
    end
    if(exist('cur_lfp_26','var'))
    per_lfp_26(isession) = average(cur_lfp_26);
    end
    if(exist('cur_lfp_34','var'))
    per_lfp_34(isession) = average(cur_lfp_34);
    end
    if(exist('cur_lfp_48','var'))
    per_lfp_48(isession) = average(cur_lfp_48);
    end
    if(exist('cur_lfp_49','var'))
    per_lfp_49(isession) = average(cur_lfp_49);
    end
    if(exist('cur_lfp_50','var'))
    per_lfp_50(isession) = average(cur_lfp_50);
    end
end
mean_pf_26=mean(pf_26_avg);
mean_pf_28=mean(pf_28_avg);
mean_pf_34=mean(pf_34_avg);
mean_pf_48=mean(pf_48_avg);
mean_pf_49=mean(pf_49_avg);
mean_pf_50=mean(pf_50_avg);
mean_lfp_28 = average(per_lfp_28);
mean_lfp_26 = average(per_lfp_26);
mean_lfp_34 = average(per_lfp_34);
mean_lfp_48 = average(per_lfp_48);
mean_lfp_49 = average(per_lfp_49);
mean_lfp_50 = average(per_lfp_50);

%function to calculate average.
function meanmatrix = average(data)
    n=length(data);
    NeuronS = zeros(101,52);
    NeuronSAvg = zeros(101,52);
    Neuront = zeros(1,101);
    for i=1:n
        [r,c]=size(data(i).S);
        if(r~=0 && c~=0)
            NeuronS(1:r,:) = data(i).S + NeuronS(1:r,:);
            NeuronSAvg(1:r,:) = 1 + NeuronSAvg(1:r,:);
            p=max(find(data(i).t>0));
            if i == 1 || ~exist('oldtlen','var')
                oldtlen = p;
                Neuront = data(i).t;
            elseif oldtlen < p
                oldtlen = p;
                Neuront = data(i).t;
            end
        end
    end
    Neuronf=data(i).f;
    AvgS=NeuronS./NeuronSAvg;
    AvgS(isnan(AvgS))=0;
    meanmatrix=struct('S',AvgS,'t',Neuront,'f',Neuronf);
end