# clean_YZ_to_4Mu_phase2_ntuplemaker

mkdir < workArea >    
cd < workArea >  
cmsrel CMSSW_10_6_27 #For UL analyses  
cd CMSSW_10_6_27/src  
git clone git@github.com:MaryHadley/clean_YZ_to_4Mu_phase2_ntuplemaker.git  
cd clean_YZ_to_4Mu_phase2_ntuplemaker  
mv * ..  
cd ..  
rm -rf clean_YZ_to_4Mu_phase2_ntuplemaker  
cmsenv

#To run Y+Z-->4 Mu code
root -l -b
.L 

