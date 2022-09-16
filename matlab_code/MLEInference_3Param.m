function [] = MLEInference_3Param(AllDatDir, RatesDir, StartChr, EndChr)
    clc;
    close all;
    if ~exist(RatesDir, 'dir')
       mkdir(RatesDir)
    end

    EnoughReads0=10; %number of reads required at t=0
    EnoughReadsLater=5; %number of reads required for later timepoints

    for chromosome = StartChr : EndChr
        inputReadDataPath = strcat(AllDatDir,'/AllDat_chr',int2str(chromosome),'.mat'); % file path of read data
        inferredRateSavingPath = strcat(RatesDir, '/Rates_chr',int2str(chromosome),'.mat');
        [fittedSites, inferedRates, inferredMethyFrac, AllDat, fitness] = kineticRatesMLE_3Param(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater);
    end
end

% clc;
% close all;
% AllDatDir = '../data/AllDat';
% RatesDir = '../data/Rates_3Param';
% StartChr = 1;
% EndChr = 1;
% if ~exist(RatesDir, 'dir')
%    mkdir(RatesDir)
% end
% 
% 
% EnoughReads0=10; %number of reads required at t=0
% EnoughReadsLater=5; %number of reads required for later timepoints
% 
% for chromosome = StartChr : EndChr
%     inputReadDataPath = strcat(AllDatDir,'/AllDat_chr',int2str(chromosome),'.mat'); % file path of read data
%     inferredRateSavingPath = strcat(RatesDir, '/Rates_chr',int2str(chromosome),'.mat');
%     [fittedSites, inferedRates, inferredMethyFrac, AllDat, Fitness] = kineticRatesMLE_3Param(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater);
% end

function [fittedSites, inferedRates, inferredMethyFrac, AllDat, Fitness] = kineticRatesMLE_3Param(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater)
        load(inputReadDataPath, 'AllDat', 'sites');
        disp(sprintf('Inferring for %s now ...', inputReadDataPath));
        %AllDat is an array of size (NSites,NTimepoints,2). The data in
        %AllDat(i,j,1) is the number of methylated reads at site i at timepoint j.
        %The data in AllDat(i,j,2) is the number of unmethylated reads at site i at
        %timepoint j. Therefore AllDat(i,j,1)+AllDat(i,j,2)=Totalreads(i,j).

        %the experimental timepoints. They are shifted by tshift, which accounts
        %for the fact that at, e.g., time=0 hrs, the captured reads started
        %replicating within the previous hour. This is done also because the model
        %assumes zero probability of methylation occuring by time=0.
        tshift=0.5;
        Times=[0, 1, 4, 16]+tshift;

        %Select only the sites which have 'enough' data--a certain number of
        %collected reads and timepoints
        Reads=sum(AllDat(:,:,1:2),3);
        NumReads0=Reads(:,1);
        NumReadsLater=sum(Reads(:,2:end),2);
        KeepSites=int64(find(NumReadsLater>=EnoughReadsLater & NumReads0>=EnoughReads0));
        numsites=int64(numel(KeepSites));

        logmink1=-2; %log10 minimum allowed fit rate
        logmaxk1=1;  %log10 maximume allowed fit rate--calculated according to ReadDepth
        logmink2=-5;
        logmaxk2=1;

        LB=[0,10^logmink1,10^logmink2];
        UB=[1,10^logmaxk1,10^logmaxk2];

        %initialize various arrays
        inferedRates=zeros(numsites,2);
        inferredMethyFrac=zeros(numsites,2);
        sites=int64(sites);
        fittedSites=int64(zeros(numsites,1));
        Fitness=zeros(numsites,1);

        tic
        SaveEvery=5E4;
        counter=0;
        x0=[1,1,0];
        t_ss=16.5; %timepoint that is considered "steady-state"

        parfor ii=1:numsites %loop over the sites
            %get the read data for this site
            Meths=AllDat(KeepSites(ii),:,1);
            UMeths=AllDat(KeepSites(ii),:,2);
            [x_fmin,f_fmin]=myfmincon(x0,Meths,UMeths,LB,UB,Times);
            
            k1=x_fmin(2);
            f=x_fmin(1);
            k2=x_fmin(3);
            
            CheckEdge=7;%10.^logmaxk1-eps
              
            if k1>=CheckEdge
                [x_CI75,f_fmin2]=myfminconCI(k1,f_fmin,f,k2,Meths,UMeths,Times);
                k1=x_CI75;
            end
            
            c1=-f*k1/(k1-k2);
            c2=f*k1/(k1-k2);
            f_ss=c1*exp(-k1*t_ss)+c2*exp(-k2*t_ss);
            Pmeth=c1*exp(-k1*Times)+c2*exp(-k2*Times);
            Pumeth=1-Pmeth;
            GetLogLikelihood=[Meths.*log(Pmeth),UMeths.*log(Pumeth)];
            Fitness(ii,:)=-sum(GetLogLikelihood,'omitnan');   
            
            inferredMethyFrac(ii,:)=[f_ss,f];
            inferedRates(ii,:)=[k1, k2];
            
            fittedSites(ii)=sites(KeepSites(ii));
        end
        save(inferredRateSavingPath,'fittedSites','inferedRates','inferredMethyFrac','Fitness')
        toc
    end

    function [x,fval]=myfmincon(x0,Meths,UMeths,LB,UB,Times)
        fun_fmin=@(x) negLL_fmin(x,Meths,UMeths,Times);
        options=optimset('Display','off');
        [x,fval,exitflag]=fmincon(fun_fmin,x0,[],[],[],[],LB,UB,[],options);
    end

    function NegLogLikelihood=negLL_fmin(ps,Meths,UMeths,Times)
        f=ps(1);
        k1=ps(2);
        k2=ps(3);
        c1=-f*k1/(k1-k2);
        c2=f*k1/(k1-k2);
        Pmeth=c1*exp(-k1*Times)+c2*exp(-k2*Times);
        Pumeth=1-Pmeth;
        GetLogLikelihood=[Meths.*log(Pmeth),UMeths.*log(Pumeth)];
        NegLogLikelihood=-sum(GetLogLikelihood,'omitnan');
    end

    function [x,fval]=myfminconCI(x0,f0,f,k2,Meths,UMeths,Times)
        fun_fmin=@(x) negLL_fminCI(x,Meths,UMeths,f0,f,k2,Times);
        options=optimset('Display','off');
        [x,fval,exitflag]=fmincon(fun_fmin,x0,[],[],[],[],[],[],[],options);
    end

    function fitness=negLL_fminCI(ps,Meths,UMeths,f0,f,k2,Times)
        PrcValIn=1.32; %75th
        k1=ps;
        c1=-f*k1/(k1-k2);
        c2=f*k1/(k1-k2);
        Pmeth=c1*exp(-k1*Times)+c2*exp(-k2*Times);
        Pumeth=1-Pmeth;
        GetLogLikelihood=[Meths.*log(Pmeth),UMeths.*log(Pumeth)];
        ProfileLogLikelihood=sum(GetLogLikelihood,'omitnan');
        newfunc=ProfileLogLikelihood+f0+PrcValIn/2;
        fitness=abs(newfunc);
    end




