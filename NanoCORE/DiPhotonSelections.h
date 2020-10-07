#ifndef DiPhotonSELECTIONS_H
#define DiPhotonSELECTIONS_H
#include "Nano.h"
#include "Base.h"


struct Photon {
    Photon(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.Photon_pt()[idx_];
        eta_ = nt.Photon_eta()[idx_];
        phi_ = nt.Photon_phi()[idx_];
        p4_ = nt.Photon_p4()[idx_];
        id_ = nt.Photon_pdgId()[idx_];
        r9_ = nt.Photon_r9()[idx_];
        chargedHadIso_ = nt.Photon_chargedHadronIso()[idx_];
        hoe_ = nt.Photon_hoe()[idx_];
        phoIso_ = nt.Photon_photonIso()[idx_];
        trkIso_ = nt.Photon_trkSumPtHollowConeDR03()[idx_];
        sieie_ = nt.Photon_sieie()[idx_];
        eveto_ = nt.Photon_electronVeto()[idx_];
        egPhoId_ = nt.Photon_mvaID()[idx_];
        //idlevel_ = whichPhotonLevel(id_, idx_);
        fixedGridRhoAll_ = nt.fixedGridRhoAll();
    }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    //int idlevel() { return idlevel_; }
    LorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    float r9() { return r9_; }
    float chargedHadIso() { return chargedHadIso_; }
    float hoe() { return hoe_; }
    float phoIso() { return phoIso_; }
    float trkIso() { return trkIso_; }
    float sieie() { return sieie_; }
    bool eveto() { return eveto_; }
    float egPhoId() { return egPhoId_; }
    float perEvtRho() { return fixedGridRhoAll_; }

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    float phi_ = 0.;
    LorentzVector p4_;
    unsigned int idx_;
    float r9_ = 0.;
    float chargedHadIso_ = 0.;
    float hoe_ = 0.;
    float phoIso_ = 0.;
    float trkIso_ = 0.;
    float sieie_ = 0.;
    bool eveto_ = 0.;
    float egPhoId_ = 0.;
    float fixedGridRhoAll_ = 0.; // this variable is the same for each event
    //int idlevel_ = SS::IDdefault;
};

vector<Photon> getPhotons();
typedef std::vector<Photon> Photons;
typedef std::pair<Photon, Photon> DiPhoton;

bool sortByPt(Photon &p1, Photon &p2)
{
        return p1.pt() > p2.pt();
        
}

vector<DiPhoton> DiPhotonPreselection(Photons &photons, bool verbose); 
bool UseLowR9Photon(Photon pho, bool isEB);
/*
SS::IDLevel whichLeptonLevel(int id, int idx);

typedef std::pair<Lepton, Lepton> Hyp;
typedef std::vector<Lepton> Leptons;

std::ostream &operator<<(std::ostream &os, Lepton &lep) {
    std::string lepstr = (abs(lep.id()) == 11) ? "Electron" : "Muon";
    return os << "<" << lepstr << " id=" << std::showpos << setw(3) << lep.id() << std::noshowpos << ", idx=" << setw(2)
              << lep.idx() << ", level=" << lep.idlevel() << ", (pT,eta)="
              << "(" << lep.pt() << "," << lep.eta() << ")>";
}
template <typename T1, typename T2> std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> &p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

vector<Lepton> getLeptons();
std::tuple<int, int, float> getJetInfo(vector<Lepton> &leps, int variation = 0);
std::pair<int, int> makesResonance(Leptons &leps, Lepton lep1, Lepton lep2, float mass, float window);
std::pair<int, Hyp> getBestHyp(vector<Lepton> &leptons, bool verbose);
bool isLeptonLevel(SS::IDLevel idlevel, int id, int idx);
void dumpLeptonProperties(Lepton lep);
*/
#endif
