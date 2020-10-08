%Note: This script needs the folders sessionInfo and MatFiles to be added to current path



clear all
clc
'This script gives spectrum and spectrogram for both spatial and feature based attention in spiral motion stimulus in 86 neurons'

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
inp_ss = input('Enter start session, please: ');
inp_es = input('Enter stop  session, please: ');
spS=0;
spF=0;
liS=0;
liF=0;
for isession = inp_ss:inp_es
    
    sessionID = ValidLFPSessions(isession);
    strcat('loading session # ', num2str(isession), '...')
    load(char(strcat(sessionID,'.mat')));
    'loaded!'
    fields=fieldnames(sessionsInfo);
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
%                     trial_data=detrend(trial_data);
%                     [P,F,Ferr]=mtspectrumc(trial_data(500:750),params);
                    trial_data=trial_data((tstart:tend));
                    %normalizing data by max
                    trial_data=trial_data/max(trial_data);
                    [S,t,f,Serr]=mtspecgramc(trial_data,movingwin,params);
                    %comment line 103 and uncomment line 100 and 106-109 for
                    %baseline normalization
%                     for i=1:length(f)
%                         baseline_p=interpl(F,P,f(i));
%                         S(:,i)=10*log10(S(:,i)/baseline_p);
%                     end
                    x = struct('S',S);
                    %unnormalized power spectrum
                    [fft1,freq,err]=mtspectrumc(trial_data,params);
                    %uncomment lines below normalization of power spectrum wrt baseline
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
            end
        end
    end
    %difference matrix for feature based
    if(exist('pf_28','var')&& exist('pf_26','var'))
        spF=spF+1;
        diff_spF(spF,1:410)=mean(pf_28)-mean(pf_26);
        per_lfp_spF(spF)=struct('S',average(cur_lfp_28).S-average(cur_lfp_26).S);
    end
    %difference matrix for spatial attention
    if(exist('pf_34','var')&& exist('pf_26','var'))
        spS=spS+1;
        diff_spS(spS,1:410)=mean(pf_34)-mean(pf_26);
        per_lfp_spS(spS)=struct('S',average(cur_lfp_34).S-average(cur_lfp_26).S);
    end
end

%spectrum for spiral motion for spatial and feature based
mean_pf_spS=mean(diff_spS);
mean_pf_spF=mean(diff_spF);

%spectrogram for spiral motion for spatial and feature based
mean_lfp_spS = average(per_lfp_spS);
mean_lfp_spF = average(per_lfp_spF);

%function to calculate average
function meanmatrix=average(data)
    [m,n]=size(data(1).S);
    S_mat=zeros(m,n);
    for i=1:length(data)
        S_mat=S_mat+data(i).S;
    end
    meanmatrix=struct('S',S_mat/length(data));
end