#include "DiPhotonSelections.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

using namespace tas;

bool UseLowR9Photon(Photon pho, bool isEB) {
    bool useThisPhoton = false;
    if (isEB) {
        if ( !(pho.sieie() < 0.015) ) return useThisPhoton;       
        if ( !(pho.trkIso() < 6.0) ) return useThisPhoton;       
        if ( !(pho.phoIso() - 0.16544*pho.perEvtRho() < 4.0) ) return useThisPhoton;       
    } else {
        if ( !(pho.sieie() < 0.035) ) return useThisPhoton;       
        if ( !(pho.trkIso() < 6.0) ) return useThisPhoton;       
        if ( !(pho.phoIso() - 0.13212*pho.perEvtRho() < 4.0) ) return useThisPhoton;       

    }
    // 0.16544 and 0.13212 are copied from flashggPreselectedDiPhotons_cfi.py
    useThisPhoton = true;
    return useThisPhoton;
}

Photons getPhotons() {
    Photons photons;
    for (unsigned int ipho = 0; ipho < nt.nPhoton(); ipho++) {
        Photon pho = Photon(ipho);
        //cout << "eta: " << pho.eta() << ", r9: " << pho.r9() << endl;
        if (pho.eveto() < 0.5) continue;
        bool pho_EB_highR9 = abs(pho.eta()) < 1.5 && pho.r9() > 0.85; 
        bool pho_EE_highR9 = abs(pho.eta()) > 1.5 && pho.r9() > 0.90; 
        bool pho_EB_lowR9 = abs(pho.eta()) < 1.5 && pho.r9() < 0.85 && pho.r9() > 0.50 && UseLowR9Photon(pho, true);
        bool pho_EE_lowR9 = abs(pho.eta()) > 1.5 && pho.r9() < 0.90 && pho.r9() > 0.80 && UseLowR9Photon(pho, false);
        if ( !(pho_EB_highR9 || pho_EE_highR9 || pho_EB_lowR9 || pho_EE_lowR9) ) continue;
        photons.push_back(pho);
    }

    sort(photons.begin(), photons.end(), sortByPt);

    //cout << "check photons pt" << endl;
    //for (unsigned int ipho = 0; ipho < nt.nPhoton(); ipho++) {
    //    cout << "pho " << ipho << ": " << photons[ipho].pt() << endl;
    //}
    return photons;
}

vector<DiPhoton> DiPhotonPreselection(Photons &photons, bool verbose=false) {
    vector<DiPhoton> diphotons; 
    for (unsigned int i1 = 0; i1 < photons.size(); i1++) {
        for (unsigned int i2 = i1+1; i2 < photons.size(); i2++) {
            Photon pho1 = photons[i1];
            Photon pho2 = photons[i2];
            DiPhoton dipho = pair<Photon, Photon>(pho1,pho2);
            if ( !(dipho.first.pt() > 35.0 && dipho.second.pt() > 25.0) ) continue;
            if ( !(dipho.first.hoe() < 0.08 && dipho.second.hoe() < 0.08) ) continue;
            if ( !(abs(dipho.first.eta()) < 2.5 && abs(dipho.second.eta()) < 2.5) ) continue;
            if ( !(abs(dipho.first.eta()) < 1.4442 || abs(dipho.first.eta()) > 1.566) ) continue;
            if ( !(abs(dipho.second.eta()) < 1.4442 || abs(dipho.second.eta()) > 1.566) ) continue;
            //if ( !((dipho.first.r9() > 0.8 && dipho.first.chargedHadIso() < 20) || dipho.first.chargedHadIso()/dipho.first.pt() < 0.3) ) continue;
            //if ( !((dipho.second.r9() > 0.8 && dipho.second.chargedHadIso() < 20) || dipho.second.chargedHadIso()/dipho.second.pt() < 0.3) ) continue;
            //if ( !(dipho.first.egPhoId() > -0.7 && dipho.second.egPhoId() > -0.7) ) continue;
            if ( !(dipho.first.r9() > 0.8 || dipho.first.chargedHadIso() < 20 || dipho.first.chargedHadIso()/dipho.first.pt() < 0.3) ) continue;
            if ( !(dipho.second.r9() > 0.8 || dipho.second.chargedHadIso() < 20 || dipho.second.chargedHadIso()/dipho.second.pt() < 0.3) ) continue;
            diphotons.push_back(dipho);
        }
    } 
    return diphotons;    
}

/*
std::tuple<int, int, float> getJetInfo(Leptons &leps, int variation) {
    int njets = 0;
    float ht = 0;
    int nbtags = 0;
    auto jetpts = Jet_pt();
    vector<float> discs = Jet_btagDeepB();
    for (unsigned int ijet = 0; ijet < jetpts.size(); ijet += 1) {
        float pt = jetpts[ijet];
        if (!(Jet_jetId()[ijet] & (1 << 1))) { continue; }
        // Clean against up to 2 fakable electrons
        if (Jet_electronIdx1()[ijet] >= 0) {
            bool skip = false;
            for (auto &lep : leps) {
                if (lep.is_el() && (Jet_electronIdx1()[ijet] == (int)(lep.idx())) && lep.idlevel() >= SS::IDfakable) {
                    skip = true;
                    break;
                }
                if (skip) { break; }
                if (Jet_electronIdx2()[ijet] >= 0) {
                    if (lep.is_el() && (Jet_electronIdx2()[ijet] == (int)(lep.idx())) && lep.idlevel() >= SS::IDfakable) {
                        skip = true;
                        break;
                    }
                    if (skip) { break; }
                }
            }
            if (skip) { continue; }
        }
        // Clean against up to 2 fakable muons
        if (Jet_muonIdx1()[ijet] >= 0) {
            bool skip = false;
            for (auto &lep : leps) {
                if (lep.is_mu() && (Jet_muonIdx1()[ijet] == (int)(lep.idx())) && lep.idlevel() >= SS::IDfakable) {
                    skip = true;
                    break;
                }
                if (skip) { break; }
                if (Jet_muonIdx2()[ijet] >= 0) {
                    if (lep.is_mu() && (Jet_muonIdx2()[ijet] == (int)(lep.idx())) && lep.idlevel() >= SS::IDfakable) {
                        skip = true;
                        break;
                    }
                    if (skip) { break; }
                }
            }
            if (skip) { continue; }
        }
        if (fabs(Jet_eta()[ijet]) > 2.4) { continue; }
        if (pt > 25. && discs[ijet] > 0.4941) { nbtags += 1; }
        if (pt < 40) { continue; }
        ht += pt;
        njets++;
    }
    return std::make_tuple(njets, nbtags, ht);
}

std::pair<int, int> makesResonance(Leptons &leps, Lepton lep1, Lepton lep2, float mass, float window) {
    // return {which lepton (1,2), and index of resonance partner}
    for (auto &lep : leps) {
        if (lep.is_el()) {
            if (!(lep1.is_el() || lep2.is_el())) { continue; }
            if ((lep.idx() == lep1.idx() && lep1.is_el()) || (lep.idx() == lep2.idx() && lep2.is_el())) { continue; }
            if (fabs(lep.eta()) > 2.4) { continue; }
            if (lep.pt() < 7) { continue; }
            if (lep.idlevel() < SS::IDveto) { continue; }
            if (lep1.is_el() && (lep1.id() * lep.id() < 0) && (fabs((lep1.p4() + lep.p4()).M() - mass) < window)) {
                return {1, lep.idx()};
            }
            if (lep2.is_el() && (lep2.id() * lep.id() < 0) && (fabs((lep2.p4() + lep.p4()).M() - mass) < window)) {
                return {2, lep.idx()};
            }
        } else {
            if (!(lep1.is_mu() || lep2.is_mu())) { continue; }
            if ((lep.idx() == lep1.idx() && lep1.is_mu()) || (lep.idx() == lep2.idx() && lep2.is_mu())) { continue; }
            if (fabs(lep.eta()) > 2.4) { continue; }
            if (lep.pt() < 5) { continue; }
            if (lep.idlevel() < SS::IDveto) { continue; }
            if (lep1.is_mu() && (lep1.id() * lep.id() < 0) && (fabs((lep1.p4() + lep.p4()).M() - mass) < window)) {
                return {1, lep.idx()};
            }
            if (lep2.is_mu() && (lep2.id() * lep.id() < 0) && (fabs((lep2.p4() + lep.p4()).M() - mass) < window)) {
                return {2, lep.idx()};
            }
        }
    }
    for (unsigned int iel = 0; iel < Electron_pt().size(); iel++) {
        if ((iel == lep1.idx() && lep1.is_el()) || (iel == lep2.idx() && lep2.is_el())) continue;
        if (fabs(Electron_eta()[iel]) > 2.4) continue;
        if (fabs(Electron_pt()[iel]) < 7) continue;
        if (!SS::electronID(iel, SS::IDveto, year())) continue;
        if (lep1.is_el() && (lep1.id() * Electron_pdgId()[iel] < 0) &&
            (fabs((lep1.p4() + Electron_p4()[iel]).M() - mass) < window)) {
            return {1, iel};
        }
        if (lep2.is_el() && (lep2.id() * Electron_pdgId()[iel] < 0) &&
            (fabs((lep2.p4() + Electron_p4()[iel]).M() - mass) < window)) {
            return {2, iel};
        }
    }
    for (unsigned int imu = 0; imu < Muon_pt().size(); imu++) {
        if ((imu == lep1.idx() && lep1.is_mu()) || (imu == lep2.idx() && lep2.is_mu())) continue;
        if (fabs(Muon_eta()[imu]) > 2.4) continue;
        if (fabs(Muon_pt()[imu]) < 5) continue;
        if (!SS::muonID(imu, SS::IDveto, year())) continue;
        if (lep1.is_mu() && (lep1.id() * Muon_pdgId()[imu] < 0) &&
            (fabs((lep1.p4() + Muon_p4()[imu]).M() - mass) < window)) {
            return {1, imu};
        }
        if (lep2.is_mu() && (lep2.id() * Muon_pdgId()[imu] < 0) &&
            (fabs((lep2.p4() + Muon_p4()[imu]).M() - mass) < window)) {
            return {2, imu};
        }
    }
    return {-1, -1};
}

std::pair<int, Hyp> getBestHyp(Leptons &leptons, bool verbose) {
    int hyp_class = -1;
    Hyp best_hyp;
    if (leptons.size() < 2) return {hyp_class, best_hyp};
    vector<Hyp> hyp1s;
    vector<Hyp> hyp2s;
    vector<Hyp> hyp3s;
    vector<Hyp> hyp4s;
    vector<Hyp> hyp6s;
    for (unsigned int i = 0; i < leptons.size(); i++) {
        for (unsigned int j = i + 1; j < leptons.size(); j++) {
            // DEBUG
            auto &lep1 = leptons[i];
            auto &lep2 = leptons[j];
            if (verbose) {
                cout << "I can read the gconf year. It's " << year() << " ok?? Get off my back." << endl;
                bool DEBUG_isss = lep1.charge() == lep2.charge();
                auto DEBUG_z_result = makesResonance(leptons, lep1, lep2, 91., 15.);
                auto DEBUG_gammastar_result = makesResonance(leptons, lep1, lep2, 0., 12.);
                bool DEBUG_extraZ = DEBUG_z_result.first >= 0;
                bool DEBUG_extraGammaStar = DEBUG_gammastar_result.first >= 0;
                cout << "hyp " << i << " leptons: " << lep2.id() << " " << lep2.pt() << " (idx: " << lep2.idx() << ") "
                     << lep1.id() << " " << lep1.pt() << " (idx: " << lep1.idx() << ")" << endl;
                cout << "   isss: " << DEBUG_isss << endl;
                cout << "   extraZ: " << DEBUG_extraZ << endl;
                cout << "   extraG: " << DEBUG_extraGammaStar << endl;
                cout << "   invt mass: " << (lep2.p4() + lep1.p4()).M() << endl;
                cout << "   passes eta: "
                     << ((lep2.is_el() ? fabs(lep2.eta()) < 2.5 : fabs(lep2.eta()) < 2.4) &&
                         (lep1.is_el() ? fabs(lep1.eta()) < 2.5 : fabs(lep1.eta()) < 2.4))
                     << " etas are " << lep2.eta() << " and " << lep1.eta() << endl;
                cout << "   passes hypsFromFirstGoodVertex: "
                     << "missing?" << endl;
                cout << "   lepton with pT " << lep2.pt() << " passes numer,denom id: " << (lep2.idlevel() == SS::IDtight)
                     << "," << (lep2.idlevel() >= SS::IDfakable) << endl;
                cout << "   lepton with pT " << lep1.pt() << " passes numer,denom id: " << (lep1.idlevel() == SS::IDtight)
                     << "," << (lep1.idlevel() >= SS::IDfakable) << endl;
                cout << "   lowMassVeto: " << ((lep2.p4() + lep1.p4()).M() < 8) << endl;
            }
            if (lep1.idlevel() < SS::IDfakable || lep2.idlevel() < SS::IDfakable) { continue; }
            bool isss = lep1.charge() == lep2.charge();
            int ntight = (lep1.idlevel() == SS::IDtight) + (lep2.idlevel() == SS::IDtight);
            int nloose = (lep1.idlevel() == SS::IDfakable) + (lep2.idlevel() == SS::IDfakable);
            if (lep1.is_el()) {
                if (lep1.pt() < 15. || fabs(lep1.eta()) > 2.5) { continue; }
            } else {
                if (lep1.pt() < 10. || fabs(lep1.eta()) > 2.4) { continue; }
            }
            if (lep2.is_el()) {
                if (lep2.pt() < 15. || fabs(lep2.eta()) > 2.5) { continue; }
            } else {
                if (lep2.pt() < 10. || fabs(lep2.eta()) > 2.4) { continue; }
            }
            // Veto SS ee or any OSSF, with mll < 12
            if ((isss && lep1.is_el() && lep2.is_el()) || (lep1.id() == -lep2.id())) {
                if ((lep1.p4() + lep2.p4()).M() < 12.) { continue; }
            }
            auto z_result = makesResonance(leptons, lep1, lep2, 91., 15.);
            auto gammastar_result = makesResonance(leptons, lep1, lep2, 0., 12.);
            bool extraZ = z_result.first >= 0;
            bool extraGammaStar = gammastar_result.first >= 0;
            int DEBUG_hyp_class = -1;
            if ((extraZ || extraGammaStar) && isss) {
                DEBUG_hyp_class = 6;
                if (lep1.pt() > lep2.pt())
                    hyp6s.push_back({lep1, lep2});
                else
                    hyp6s.push_back({lep2, lep1});
            } else if (ntight == 2 && isss) {
                DEBUG_hyp_class = 3;
                if (lep1.pt() > lep2.pt())
                    hyp3s.push_back({lep1, lep2});
                else
                    hyp3s.push_back({lep2, lep1});
            } else if (nloose == 2 && isss) {
                DEBUG_hyp_class = 1;
                if (lep1.pt() > lep2.pt())
                    hyp1s.push_back({lep1, lep2});
                else
                    hyp1s.push_back({lep2, lep1});
            } else if (ntight == 1 && nloose == 1 && isss) {
                DEBUG_hyp_class = 2;
                if (lep1.pt() > lep2.pt())
                    hyp2s.push_back({lep1, lep2});
                else
                    hyp2s.push_back({lep2, lep1});
            } else if (ntight == 2 && !isss) {
                DEBUG_hyp_class = 4;
                if (lep1.pt() > lep2.pt())
                    hyp4s.push_back({lep1, lep2});
                else
                    hyp4s.push_back({lep2, lep1});
            }
            if (verbose) cout << "hyp #" << i << " hyp_class: " << DEBUG_hyp_class << endl;
        }
    }
    vector<Hyp> hyps;
    // Priority order - 3,6,2,1,4
    if (hyp3s.size() > 0) {
        hyps = hyp3s;
        hyp_class = 3;
    } else if (hyp6s.size() > 0) {
        hyps = hyp6s;
        hyp_class = 6;
    } else if (hyp2s.size() > 0) {
        hyps = hyp2s;
        hyp_class = 2;
    } else if (hyp1s.size() > 0) {
        hyps = hyp1s;
        hyp_class = 1;
    } else if (hyp4s.size() > 0) {
        hyps = hyp4s;
        hyp_class = 4;
    }
    if ((hyp_class <= 0) || (hyps.size() < 1)) return {hyp_class, best_hyp};
    if (hyps.size() == 1) {
        best_hyp = hyps[0];
    } else if (hyps.size() > 1) {
        best_hyp = hyps[0];
        for (unsigned int i = 1; i < hyps.size(); i++) {
            Hyp hyp = hyps[i];
            if (hyp.first.is_mu() + hyp.second.is_mu() > best_hyp.first.is_mu() + best_hyp.second.is_mu()) {
                best_hyp = hyp;
            } else if (hyp.first.is_mu() + hyp.second.is_mu() == best_hyp.first.is_mu() + best_hyp.second.is_mu()) {
                if (hyp.first.pt() + hyp.second.pt() > best_hyp.first.pt() + best_hyp.second.pt()) { best_hyp = hyp; }
            }
        }
    }
    return {hyp_class, best_hyp};
}

void dumpLeptonProperties(Lepton lep) {
    std::cout << lep << std::endl;
    int i = lep.idx();
    if (lep.is_el()) {
        std::cout << "  -- etaSC: " << Electron_eta()[i] + Electron_deltaEtaSC()[i] << std::endl;
        std::cout << "  -- mva: " << Electron_mvaFall17V1noIso()[i] << std::endl;
        std::cout << "  -- lostHits: " << Electron_lostHits()[i] << std::endl;
        std::cout << "  -- miniRelIso: " << Electron_miniPFRelIso_all()[i] << std::endl;
        std::cout << "  -- ptRatio: " << 1/(Electron_jetRelIso()[i] + 1) << std::endl;
        std::cout << "  -- ptRel: " << Electron_jetPtRelv2()[i] << std::endl;
    } else {
        std::cout << "  -- miniRelIso: " << Muon_miniPFRelIso_all()[i] << std::endl;
        std::cout << "  -- ptRatio: " << 1/(Muon_jetRelIso()[i] + 1) << std::endl;
        std::cout << "  -- ptRel: " << Muon_jetPtRelv2()[i] << std::endl;
    }
}

bool isLeptonLevel(SS::IDLevel idlevel, int id, int idx) {
    if (abs(id) == 11) {
        return SS::electronID(idx, idlevel, year());
    } else if (abs(id) == 13) {
        return SS::muonID(idx, idlevel, year());
    } else {
        return false;
    }
}

SS::IDLevel whichLeptonLevel(int id, int idx) {
    if (isLeptonLevel(SS::IDveto, id, idx)) {
        if (isLeptonLevel(SS::IDfakable, id, idx)) {
            if (isLeptonLevel(SS::IDtight, id, idx)) {
                return SS::IDtight;
            } else {
                return SS::IDfakable;
            }
        } else {
            return SS::IDveto;
        }
    } else {
        return SS::IDdefault;
    }
}
*/
