%Note: This script needs the folders sessionInfo and MatFiles to be added to current path

clear all
clc

'This script gives spectrum and spectrogram for both spatial and feature based attention in both linear and spiral motion stimulus in 34 neurons'

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
N=800;
nyquist=srate/2;
frequencies = linspace(0,nyquist,floor(N/2)+1);
inp_ss = input('Enter start session, please: ');
inp_es = input('Enter stop  session, please: ');
spS=0;
spF=0;
liS=0;
liF=0;
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
    sessionInfo = sessionsInfo(session).sessionInfo; %Information on the timing of each trial
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
%                     [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    %normalization by max
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %comment line 117 and uncomment line 114 and 120-123 for baseline
                    %normalization
%                     for i=1:length(f)
%                         baseline_p=interpl(F,P,f(i));;
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    x = struct('S',S);
                    [fft1,freq,err]=mtspectrumc(trial_data,params);
                    %uncomment lines below to normalize spectrum wrt
                    %baseline
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
%                     [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
%                   comment line 153 and uncomment line 151 and 156-159 for baseline normalization of spectrogram 
%                     for i=1:length(f)
%                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    x = struct('S',S);
                    [fft1,freq,err]=mtspectrumc(trial_data(1:1200),params);
                    %uncomment below lines to normalize spectrum wrt
                    %baseline
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
                    %uncomment line 187 for baseline normalization
                    %[P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    
                    trial_data=trial_data((tstart:tend));
                    %normalization by max. comment line 191 for baseline
                    %normalization
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %uncomment lines below for baseline normalization of
                    %spectrogram
%                     for i=1:length(f)
%                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    x = struct('S',S);
                    [fft1,freq,err]=mtspectrumc(trial_data(1:1200),params);
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
    %difference matrix for spatial attention in LMS
    if(exist('pf_50','var')&& exist('pf_48','var'))
        liS=liS+1;
        x=size(pf_50);
        if(x(1)==1)
            diff_liS(liS,1:410)=pf_50-mean(pf_48);
        else
            diff_liS(liS,1:410)=mean(pf_50)-mean(pf_48);
        end
        per_lfp_liS(liS)=struct('S',average(cur_lfp_50).S-average(cur_lfp_48).S);
    end
    %difference matrix for feature based attention in LMS
    if(exist('pf_49','var')&& exist('pf_48','var'))
        liF=liF+1;
        diff_liF(liF,1:410)=mean(pf_48)-mean(pf_49);
        per_lfp_liF(liF)=struct('S',average(cur_lfp_48).S-average(cur_lfp_49).S);
    end
    %difference matrix for feature bases attention in SMS
    if(exist('pf_28','var')&& exist('pf_26','var'))
        spF=spF+1;
        diff_spF(spF,1:410)=mean(pf_28)-mean(pf_26);
        per_lfp_spF(spF)=struct('S',average(cur_lfp_28).S-average(cur_lfp_26).S);
    end
    %difference matrix for spatial attention in SMS
    if(exist('pf_34','var')&& exist('pf_26','var'))
        spS=spS+1;
        diff_spS(spS,1:410)=mean(pf_34)-mean(pf_26);
        per_lfp_spS(spS)=struct('S',average(cur_lfp_34).S-average(cur_lfp_26).S);
    end
end
mean_pf_spS=mean(diff_spS);
mean_pf_spF=mean(diff_spF);
mean_pf_liS=mean(diff_liS);
mean_pf_liF=mean(diff_liF);
mean_lfp_spS = average(per_lfp_spS);
mean_lfp_spF = average(per_lfp_spF);
mean_lfp_liF = average(per_lfp_liF);
mean_lfp_liS = average(per_lfp_liS);
function meanmatrix=average(data)
    [m,n]=size(data(1).S);
    S_mat=zeros(m,n);
    for i=1:length(data)
        S_mat=S_mat+data(i).S;
    end
    meanmatrix=struct('S',S_mat/length(data));
end