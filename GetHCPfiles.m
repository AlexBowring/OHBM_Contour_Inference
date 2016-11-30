OutBase='/home/essicd/storage/data/ContourInf/HCP80'
HCP='/home/essicd/storage/data/HCP/Unrelated_80/RESULTS'

Dirs=dir([HCP '/stats_*']);

for i=3:length(Dirs)
  Src=Dirs(i);
  [~,OutNm]=fileparts(Src.name);
  Out=fullfile(OutBase,OutNm);
  mkdir(Out);
  String=[HCP '/' Src.name '/'];
  BrainBootstrappingFSL(String, Out);
end
