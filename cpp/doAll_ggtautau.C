{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_ggtautau.C+");
    TChain *ch = new TChain("Events");

    //ch->Add("root://redirector.t2.ucsd.edu///store/user/namin/nanoaod/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/8E0C8306-DC0D-0548-BA7C-D0698140DF28.root");
    //ch->Add("root://redirector.t2.ucsd.edu///store/user/hmei/nanoaod_runII/HHggtautau/HHggtautau_Era2018_private_v2_20201005/test_nanoaod_1.root");
    ch->Add("root:///home/users/hmei/myWorkspace/privateMC_gen/CMSSW_10_2_22/src/HIG-RunIIAutumn18NanoAODv7-01134.root");

    ScanChain(ch);

}
