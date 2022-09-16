function [] = ExtractReversibleSites(AllDatDir, Rates2ParamDir, Rates3ParamDir, ReversibleSitesSaveDir, StartChr, EndChr)
    clc;
    close all;
%     AllDatDir = '../data/AllDat';
%     Rates2ParamDir = '../data/Rates_2Param';
%     Rates3ParamDir = '../data/Rates_3Param';
%     ReversibleSitesSaveDir = '../data/ReversibleSites';
%     StartChr = 1;
%     EndChr = 22;
    if ~exist(ReversibleSitesSaveDir, 'dir')
       mkdir(ReversibleSitesSaveDir)
    end
    for chromosome = StartChr : EndChr
        inputReadDataPath = strcat(AllDatDir,'/AllDat_chr',int2str(chromosome),'.mat'); % file path of read data
        inferred2ParamRateSavingPath = strcat(Rates2ParamDir, '/Rates_chr',int2str(chromosome),'.mat');
        inferred3ParamRateSavingPath = strcat(Rates3ParamDir, '/Rates_chr',int2str(chromosome),'.mat');
        reversibleSitesSavingPath = strcat(ReversibleSitesSaveDir, '/ReversibleSites_chr',int2str(chromosome),'.mat');
        reversibleSites = Compare2_and_3Param(inputReadDataPath, inferred2ParamRateSavingPath, inferred3ParamRateSavingPath, reversibleSitesSavingPath, chromosome);
    end
end



function reversibleSites = Compare2_and_3Param(inputReadDataPath, inferred2ParamRateSavingPath, inferred3ParamRateSavingPath, reversibleSitesSavingPath, chromosome)
	load(inputReadDataPath,'AllDat','sites');

	Reads=sum(AllDat(:,:,1:2),3);
	NumReads0=Reads(:,1);
	NumReadsLater=sum(Reads(:,2:end),2);
	NumReadsAll=NumReads0 + NumReadsLater;

	%get the raw avg. percent methylation at each timepoint
	Tots=sum(AllDat,3);
	PercMethSites=AllDat(:,:,1)./Tots;

	load(inferred3ParamRateSavingPath);
	%load in the 3-parameter fits
	P3Lam=inferedRates;
	P3Frac=inferredMethyFrac;
	P3Fitness=Fitness;
	P3Sites=fittedSites;
	clear inferedRates inferredMethyFrac Fitness fittedSites

	%identify how many reads were measured for the sites that were kept
	Lia=ismember(sites,P3Sites);
	siteinds=find(Lia); %these are the index locations in AllDat file of analyzed sites
	NumReads=NumReadsAll(siteinds);

	%load in the 2-parameter fits
	load(inferred2ParamRateSavingPath);
	P2Lam=inferedRates;
	P2Frac=inferredMethyFrac;
	P2Fitness=-Fitness;
	P2Sites=fittedSites; %should be the same as P3 sites
	clear inferedRates inferredMethyFrac Fitness fittedSites

	%compute the information criteria (AIC, BIC)
	%AIC metric is 2k-2ln(L). k=#params, L=optimized Likelihood. Our fitnesss =
	%negative log likelihood
	%AIC2=2*2+2*P2Fitness;
	%AIC3=2*3+2*P3Fitness;

	%BIC metric is ln(n)k-2ln(L). n=#data pts (NumReads)
	BIC2=log(NumReads)*2+2*P2Fitness;
	BIC3=log(NumReads)*3+2*P3Fitness;

	%identify "reversible sites" as those for which BIC3 (3 param) is sig. less
	%than BIC2 (2 param)
	DiffBIC=BIC3-BIC2;
	%find sites for which 3 param is lower
	%JustLess=find(DiffBIC<0);
	LessBy2=find(DiffBIC<-2); %Significant difference if less by 2 or more

	%Now keep just the parameters and site indices for the identified
	%"reversibel sites"
	reversibleRates=P3Lam(LessBy2,:);
	reversibleFracs=P3Frac(LessBy2,:);
	reversibleSites=P3Sites(LessBy2);
	PctReversible=numel(reversibleSites)/numel(P2Sites);
	disp(sprintf('# of reversible sites for chr%d: %d/%d (%.1f%%)', chromosome, numel(reversibleSites), numel(P2Sites), PctReversible*100.0));
	%save to file
	save(reversibleSitesSavingPath,'reversibleRates', 'reversibleFracs', 'reversibleSites');
	clear AllDat sites Tots PercMethSites Lia siteinds Reads NumReads NumReads0 NumReadsLater NumReadsAll reversibleRates reversibleFracs P3Lam P3Frac P3Sites P3Fitness P2Lam P2Frac P2Sites P2Fitness BIC2 BIC3 LessBy2 DiffBIC PctReversible
end