#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

#include <vector>
#include <iostream>

using namespace std;

enum LIMITS {
    MAX_TRACKS = 30000,
};

bool is_in_fv(const TLorentzVector &pos, const float margin = 20)
{
    if (pos.X() > 3 + margin && pos.X() < 253 - margin && pos.Y() > -113 + margin && pos.Y() < 113 - margin &&
        pos.Z() > 3 + margin && pos.Z() < 1034 - margin) {
        return true;
    }

    return false;
}

bool is_in_fv(const float x, const float y, const float z, const float margin = 20)
{
    TLorentzVector pos(x, y, z, 0);
    return is_in_fv(pos, margin);
}

bool is_good_proton(const TLorentzVector &start, const TLorentzVector &end, const float proton_daugher_ke,
                    const float max_proton_daugher_ke = 0)
{
    if (is_in_fv(start, 0) && is_in_fv(end, 0) && proton_daugher_ke <= max_proton_daugher_ke) {
        return true;
    }
    return false;
}

TH1 *resolution_hist(TH2 *in)
{
    TString basename(in->GetName());
    auto current_dir = in->GetDirectory();
    in->FitSlicesY();
    auto mean = (TH1 *) current_dir->Get(basename + "_1");
    auto sigma = (TH1 *) current_dir->Get(basename + "_2");
    auto out = (TH1 *) mean->Clone(basename + "res");
    for (int ix = 1; ix < out->GetNbinsX(); ++ix) {
        out->SetBinError(ix, sigma->GetBinContent(ix));
    }
    return out;
}

bool particle_match(const TLorentzVector &pa, const TLorentzVector &pb, const TLorentzVector &ma, const TLorentzVector &mb)
{
    // pos + direction
    // if (fabs(pa.X() - pb.X() + 0.1) > 0.5) return false;
    // if (fabs(pa.Y() - pb.Y()) > 1.) return false;
    // if (fabs(pa.Z() - pb.Z()) > 1.) return false;
    // if (fabs(ma.Theta() - mb.Theta()) > 0.5) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 0.5) return false;

    // dist
    auto dist = (pa - pb).Vect();
    if (dist.Mag() > 1.) return false;

    // if (fabs(ma.Theta() - mb.Theta()) > 1.0) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 1.0) return false;

    return true;
}

class PF {
   public:
    PF(float inu_energy, int iNtrack, int *iid, int *ipdg, float *istartMomentum, float *istartXYZT, float *iendXYZT,
       int *imother, std::vector<std::vector<int> > *idaughters)
    {
        nu_energy = inu_energy;
        Ntrack = iNtrack;
        id = iid;
        pdg = ipdg;
        startMomentum = istartMomentum;
        startXYZT = istartXYZT;
        endXYZT = iendXYZT;
        mother = imother;
        daughters = idaughters;
    }
    void print()
    {
        cout << Ntrack << endl;
        for (int i = 0; i < Ntrack; ++i) {
            cout
            << id[i] << ", "
            << pdg[i] << ", "
            << mother[i] << ", "
            << daughters->at(i).size() << ", "
            << startMomentum[i * 4 + 3] << endl;
        }
    }
    string dotify() {
        // std::map<int,char*> m_color;
        // m_color[2212] = "black";
        // m_color[-2212] = "black";
        // m_color[211] = "blue";
        // m_color[-211] = "blue";
        // m_color[13] = "chartreuse3";
        // m_color[-13] = "chartreuse3";
        // m_color[11] = "cyan2";
        // m_color[-11] = "cyan2";
        // m_color[22] = "cyan3";
        stringstream ss;
        ss << "digraph D {\n";
        ss << 0 << " [label=\"" << "neutrino" << "\\l" << nu_energy << "\"]\n";
        for (int i = 0; i < Ntrack; ++i) {
            float px = startMomentum[i * 4 + 0];
            float py = startMomentum[i * 4 + 1];
            float pz = startMomentum[i * 4 + 2];
            float e = startMomentum[i * 4 + 3];
            TLorentzVector tlv(px,py,pz,e);
            auto ke = tlv.E() - tlv.M();
            string color = "black";
            stringstream fv;
            if (!is_in_fv(startXYZT[i * 4 + 0], startXYZT[i * 4 + 1], startXYZT[i * 4 + 2], 0)) {
                fv << " s.o.f ";
            }
            if (!is_in_fv(endXYZT[i * 4 + 0], endXYZT[i * 4 + 1], endXYZT[i * 4 + 2], 0)) {
                // fv << endXYZT[i * 4 + 0] << ", " << endXYZT[i * 4 + 1] << ", " << endXYZT[i * 4 + 2];
                fv << " e.o.f ";
            }
            if(ke<0.01) continue;
            ss << id[i] << " [label=\"" << pdg[i] << "\\l" << ke << "\\l" << fv.str() << "\"";
            // ss << "color=" << color;
            ss << "]\n";
            ss << mother[i] << "->" << id[i] << endl;
        }
        ss << "}\n";
        return ss.str();
    }
    float nu_energy;
    int Ntrack;
    int *id;
    int *pdg;
    float *startMomentum;
    float *startXYZT;
    float *endXYZT;
    int *mother;
    std::vector<std::vector<int> > *daughters;
};

void proton_eff(
    // const char *input = "data/nue_overlay_run2.root"
    const char *input = "data/nu_overlay_run2.root"
    )
{
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    const int target_pdg = 2212;  // 11, 13, 211, 2212
    const float min_ke_multiplicity = 0.035; // GeV

    auto *tf = TFile::Open(input, "read");
    auto *dir = (TDirectoryFile *) tf->Get("wcpselection");

    // old
    // auto *T_PFDump = (TTree *) dir->Get("T_PFDump");
    // auto *T_PFeval = (TTree *) dir->Get("T_PFeval");
    // T_PFDump->AddFriend(T_PFeval);

    // new
    auto *T_PFDump = (TTree *) dir->Get("T_PFeval");

    auto *T_KINEvars = (TTree *) dir->Get("T_KINEvars");
    T_PFDump->AddFriend(T_KINEvars);
    auto *T_eval = (TTree *) dir->Get("T_eval");
    T_PFDump->AddFriend(T_eval);
    auto *T_BDTvars = (TTree *) dir->Get("T_BDTvars");
    T_PFDump->AddFriend(T_BDTvars);

    int run;
    int subrun;
    int event;

    int truth_Ntrack;
    int truth_id[MAX_TRACKS];
    int truth_pdg[MAX_TRACKS];
    int truth_process[MAX_TRACKS];
    int truth_mother[MAX_TRACKS];
    float truth_startXYZT[MAX_TRACKS][4];
    float truth_endXYZT[MAX_TRACKS][4];
    float truth_startMomentum[MAX_TRACKS][4];
    float truth_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *truth_daughters = 0;

    int reco_Ntrack;
    int reco_id[MAX_TRACKS];
    int reco_pdg[MAX_TRACKS];
    int reco_process[MAX_TRACKS];
    int reco_mother[MAX_TRACKS];
    float reco_startXYZT[MAX_TRACKS][4];
    float reco_endXYZT[MAX_TRACKS][4];
    float reco_startMomentum[MAX_TRACKS][4];
    float reco_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *reco_daughters = 0;

    bool truth_isCC;
    float truth_nu_momentum[4];
    float truth_corr_nuvtxX;
    float truth_corr_nuvtxY;
    float truth_corr_nuvtxZ;
    float reco_nuvtxX;
    float reco_nuvtxY;
    float reco_nuvtxZ;

    float numu_cc_flag;
    float nue_score;
    float numu_score;

    float stm_clusterlength;
    bool truth_isFC;
    bool match_isFC;

    float kine_reco_Enu;
    std::vector<float> *kine_energy_particle = 0;
    std::vector<int> *kine_energy_info = 0;
    std::vector<int> *kine_particle_type = 0;

    T_PFDump->SetBranchAddress("run", &run);
    T_PFDump->SetBranchAddress("subrun", &subrun);
    T_PFDump->SetBranchAddress("event", &event);

    T_PFDump->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_PFDump->SetBranchAddress("truth_id", &truth_id);
    T_PFDump->SetBranchAddress("truth_pdg", &truth_pdg);
    T_PFDump->SetBranchAddress("truth_mother", &truth_mother);
    T_PFDump->SetBranchAddress("truth_startXYZT", &truth_startXYZT);
    T_PFDump->SetBranchAddress("truth_endXYZT", &truth_endXYZT);
    T_PFDump->SetBranchAddress("truth_startMomentum", &truth_startMomentum);
    T_PFDump->SetBranchAddress("truth_daughters", &truth_daughters);

    T_PFDump->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_PFDump->SetBranchAddress("reco_id", &reco_id);
    T_PFDump->SetBranchAddress("reco_pdg", &reco_pdg);
    T_PFDump->SetBranchAddress("reco_mother", &reco_mother);
    T_PFDump->SetBranchAddress("reco_startXYZT", &reco_startXYZT);
    T_PFDump->SetBranchAddress("reco_endXYZT", &reco_endXYZT);
    T_PFDump->SetBranchAddress("reco_startMomentum", &reco_startMomentum);
    T_PFDump->SetBranchAddress("reco_daughters", &reco_daughters);

    T_PFDump->SetBranchAddress("truth_isCC", &truth_isCC);
    T_PFDump->SetBranchAddress("truth_nu_momentum", &truth_nu_momentum);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxX", &truth_corr_nuvtxX);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxY", &truth_corr_nuvtxY);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxZ", &truth_corr_nuvtxZ);
    T_PFDump->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);
    T_PFDump->SetBranchAddress("reco_nuvtxY", &reco_nuvtxY);
    T_PFDump->SetBranchAddress("reco_nuvtxZ", &reco_nuvtxZ);

    //  T_BDTvars
    T_PFDump->SetBranchAddress("numu_cc_flag", &numu_cc_flag);
    T_PFDump->SetBranchAddress("nue_score", &nue_score);
    T_PFDump->SetBranchAddress("numu_score", &numu_score);

    // T_eval
    T_PFDump->SetBranchAddress("stm_clusterlength", &stm_clusterlength);
    T_PFDump->SetBranchAddress("truth_isFC", &truth_isFC);
    T_PFDump->SetBranchAddress("match_isFC", &match_isFC);

    T_PFDump->SetBranchAddress("kine_energy_particle", &kine_energy_particle);
    T_PFDump->SetBranchAddress("kine_energy_info", &kine_energy_info);
    T_PFDump->SetBranchAddress("kine_particle_type", &kine_particle_type);
    T_PFDump->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu);

    TH1F *h_ndaughter_truth = new TH1F("h_ndaughter_truth", "h_ndaughter_truth", 100, -0.5, 99.5);
    TH1F *h_ndaughter_reco = new TH1F("h_ndaughter_reco", "h_ndaughter_reco", 100, -0.5, 99.5);

    // TH1F *h_reco_m_truth = new TH1F("h_reco_m_truth", "h_reco_m_truth", 10,-0.5,9.5); // multiplicity
    TH1F *h_reco_m_truth = new TH1F("h_reco_m_truth", "h_reco_m_truth", 100,-3,3); // multiplicity

    TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 100, 0, 2, 100, 0, 2); // energy
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 2, 40, 0, 2); // energy bias
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 10, 0, 3, 10, -0.5, 9.5); // multiplicity
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 100, -TMath::Pi(), TMath::Pi()); //
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 10, 0, 3, 100, -3, 3);  // angle
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 100, -1, 1);  // tmp

    TH1F *h_truth_e_all = new TH1F("h_truth_e_all", "h_truth_e_all", 60, 0, 3);
    TH1F *h_truth_e_match = new TH1F("h_truth_e_match", "h_truth_e_match", 60, 0, 3);
    // TH1F *h_truth_e_all = new TH1F("h_truth_e_all", "h_truth_e_all", 10, -0.5, 9.5);
    // TH1F *h_truth_e_match = new TH1F("h_truth_e_match", "h_truth_e_match", 10, -0.5, 9.5);

    TH2F *h_dtheta_theta =
        new TH2F("h_dtheta_theta", "h_dtheta_theta", 200, -TMath::Pi(), TMath::Pi(), 200, -TMath::Pi(), TMath::Pi());

    int counter_all = 0;
    int counter_pass = 0;
    // for (int ientry = 0; ientry < T_PFDump->GetEntries(); ++ientry) {
    for (int ientry = 0; ientry < 10000; ++ientry) {
        T_PFDump->GetEntry(ientry);
        if (ientry % (T_PFDump->GetEntries()/10) == 0) cout << "processing: " << ientry*100 / T_PFDump->GetEntries() << "%" << endl;
        // if (ientry != 30779) continue;
        if (numu_cc_flag < 0 || stm_clusterlength < 15) continue;  // generic nu selection
        if (numu_score<0.9) continue; // BDT numu
        // if (nue_score<7) continue; // BDT nue
        if (!truth_isCC) continue;
        if (match_isFC!=true) continue; // FV cut
        // if (truth_isFC!=true) continue; // FV cut
        // if (!is_in_fv(truth_corr_nuvtxX, truth_corr_nuvtxY, truth_corr_nuvtxZ)) continue; // alternative FV
        // if (truth_nu_momentum[3]<0.3 || truth_nu_momentum[3]>1.5) continue;


        TVector3 truth_nuvtx(truth_corr_nuvtxX, truth_corr_nuvtxY, truth_corr_nuvtxZ);
        TVector3 reco_nuvtx(reco_nuvtxX, reco_nuvtxY, reco_nuvtxZ);
        if ((reco_nuvtx - truth_nuvtx).Mag() > 1.0) continue;

        TLorentzVector target_pos_start_truth;
        TLorentzVector target_pos_end_truth;
        TLorentzVector target_mom_start_truth;
        int target_index_truth = -1;
        float current_max_energy_truth = -1;
        // find the leading primary target particle
        for (int itruth = 0; itruth < truth_Ntrack; ++itruth) {
            h_ndaughter_truth->Fill(truth_daughters->at(itruth).size());
            if (truth_mother[itruth] != 0) continue;
            if (abs(truth_pdg[itruth]) != 2212) continue;
            if(truth_startMomentum[itruth][3]<current_max_energy_truth) continue;
            target_index_truth = itruth;
            current_max_energy_truth = truth_startMomentum[itruth][3];
            target_pos_start_truth.SetXYZT(truth_startXYZT[itruth][0], truth_startXYZT[itruth][1],
                                           truth_startXYZT[itruth][2], truth_startXYZT[itruth][3]);
            target_pos_end_truth.SetXYZT(truth_endXYZT[itruth][0], truth_endXYZT[itruth][1], truth_endXYZT[itruth][2],
                                         truth_endXYZT[itruth][3]);
            target_mom_start_truth.SetXYZT(truth_startMomentum[itruth][0], truth_startMomentum[itruth][1],
                                           truth_startMomentum[itruth][2], truth_startMomentum[itruth][3]);
        }

        // leading reco target
        TLorentzVector pos_reco;
        TLorentzVector mom_reco;
        float current_max_energy_reco = -1;
        for (int ireco = 0; ireco < reco_Ntrack; ++ireco) {
            h_ndaughter_reco->Fill(reco_daughters->at(ireco).size());
            if (reco_mother[ireco] != 0) continue;
            if (abs(reco_pdg[ireco]) != 2212) continue;
            if (reco_startMomentum[ireco][3]<current_max_energy_reco) continue;
            current_max_energy_reco = reco_startMomentum[ireco][3];
            pos_reco.SetXYZT(reco_startXYZT[ireco][0], reco_startXYZT[ireco][1], reco_startXYZT[ireco][2],
                             reco_startXYZT[ireco][3]);
            mom_reco.SetXYZT(reco_startMomentum[ireco][0], reco_startMomentum[ireco][1], reco_startMomentum[ireco][2],
                             reco_startMomentum[ireco][3]);
        }

        if (current_max_energy_truth<0) {
            continue;
        }

        bool simple_proton = true;
        if (!is_in_fv(truth_endXYZT[target_index_truth][0], truth_endXYZT[target_index_truth][1],
                        truth_endXYZT[target_index_truth][2], 0)) {
            simple_proton = false;
        }
        for (int itruth = 0; itruth < truth_Ntrack; ++itruth) {
            if (truth_mother[itruth] != truth_id[target_index_truth]) continue;
            TLorentzVector tlv;
            tlv.SetXYZT(truth_startMomentum[itruth][0], truth_startMomentum[itruth][1],
                                           truth_startMomentum[itruth][2], truth_startMomentum[itruth][3]);
            auto ke = tlv.E() - tlv.M();
            if (ke > 0.01) simple_proton = false;
        }
        if (!simple_proton) {
            continue;
        }

        ++counter_all;
        const float target_KE = target_mom_start_truth.E()-target_mom_start_truth.M();
        h_truth_e_all->Fill(target_KE);

        if (true && target_KE > 0.1) {
            char buff[100];
            snprintf(buff, sizeof(buff), "%1.2f", target_KE);
            stringstream ss;
            if (current_max_energy_reco < 0) {
                ss << "fail_";
            } else {
                ss << "match_";
            }
            ss << buff << "_" << ientry << "_" << run << "_" << subrun << "_" << event;
            ofstream ftruth(ss.str() + "_truth.dot");
            PF pf_truth(truth_nu_momentum[3], truth_Ntrack, truth_id, truth_pdg, (float *) truth_startMomentum, (float *) truth_startXYZT, (float *) truth_endXYZT,
                        truth_mother, truth_daughters);
            ftruth << pf_truth.dotify();
            ftruth.close();
            ofstream freco(ss.str() + "_reco.dot");
            PF pf_reco(kine_reco_Enu / 1000., reco_Ntrack, reco_id, reco_pdg, (float *) reco_startMomentum, (float *) reco_startXYZT,(float *) reco_endXYZT, reco_mother,
                       reco_daughters);
            freco << pf_reco.dotify();
            freco.close();
        }

        if (current_max_energy_reco < 0) {
            continue;
        }
        ++counter_pass;
        h_truth_e_match->Fill(target_KE);
        h_reco_v_truth->Fill(target_KE,mom_reco.E()-mom_reco.M());
    }
    cout << "target: truth " << counter_all << "; reco " << counter_pass
         << "; ratio: " << 100. * counter_pass / counter_all << "%" << endl;

    const int LINE_WIDTH = 2;

    TCanvas *c0 = new TCanvas("c0", "c0");
    h_truth_e_all->SetLineColor(kBlack);
    h_truth_e_all->SetLineWidth(LINE_WIDTH);
    h_truth_e_all->SetTitle(";E^{truth} [GeV]");
    h_truth_e_all->Draw();
    h_truth_e_match->SetLineColor(kRed);
    h_truth_e_match->SetLineWidth(LINE_WIDTH);
    h_truth_e_match->Draw("same");

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    auto pEff = new TEfficiency(*h_truth_e_match, *h_truth_e_all);
    // pEff->GetXaxis()->SetRangeUser(0,1);
    pEff->SetMarkerStyle(20);
    pEff->SetTitle(";E^{truth} [GeV];Efficiency");
    // pEff->SetTitle(";True # proton (>35 MeV);Efficiency");
    pEff->Draw("ap");
    // auto h_truth_e_ratio = (TH1F *) h_truth_e_match->Clone("h_truth_e_ratio");
    // h_truth_e_ratio->Divide(h_truth_e_all);
    // h_truth_e_ratio->SetLineColor(kRed);
    // h_truth_e_ratio->SetLineWidth(LINE_WIDTH);
    // h_truth_e_ratio->SetTitle(";E^{truth} [GeV];Efficiency");
    // h_truth_e_ratio->SetStats(0);
    // h_truth_e_ratio->Draw("");

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->SetGrid();
    c2->SetLogz();
    // h_reco_v_truth->SetTitle(";Z^{truth} [cm];Z^{reco}-Z^{truth} [cm]");
    // h_reco_v_truth->SetTitle(";#phi^{truth} [rad];#phi^{reco}-#phi^{truth} [rad]");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];reco/true");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];multiplicity");
    h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco} [GeV]");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];#Delta #theta");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];#Delta #phi");
    // h_reco_v_truth->SetStats(0);
    h_reco_v_truth->SetMinimum(1);
    h_reco_v_truth->Draw("colz");

    // resolution_hist
    // auto h_reco_v_truth_1 = resolution_hist(h_reco_v_truth);
    // h_reco_v_truth_1->Draw("e,same");
    // h_reco_v_truth->QuantilesX()->Draw("same");
    // h_reco_v_truth->ProfileX()->Draw("same");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->SetGrid();
    // h_reco_m_truth->SetTitle(";|Vtx^{reco}-Vtx^{truth}| [cm]");
    // h_reco_m_truth->SetTitle(";X^{reco}-X^{truth} [cm]");
    // h_reco_m_truth->SetTitle(";#phi^{reco}-#phi^{truth} [rad]");
    h_reco_m_truth->SetTitle(";E^{reco}-E^{truth} [GeV]");
    // h_reco_m_truth->SetTitle(";# proton");
    // h_reco_m_truth->SetStats(0);
    h_reco_m_truth->Draw();

    // TCanvas *c4 = new TCanvas("c4", "c4");
    // h_ndaughter_truth->Draw();
    // TCanvas *c5 = new TCanvas("c5", "c5");
    // h_ndaughter_reco->Draw();

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->SetGrid();
    h_dtheta_theta->SetTitle("; #theta_{YZ}; #theta_{XY}");
    h_dtheta_theta->Draw("colz");
}
