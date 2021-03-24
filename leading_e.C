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

bool fv_shwr_vtx(const TLorentzVector &pos)
{
    const float margin = 20;
    if (pos.X() > 3 + margin && pos.X() < 253 - margin && pos.Y() > -113 + margin && pos.Y() < 113 - margin &&
        pos.Z() > 3 + margin && pos.Z() < 1034 - margin) {
        return true;
    }

    return false;
}

TH1 *resolution_hist(TH2 *in)
{
    TString basename(in->GetName());
    auto current_dir = in->GetDirectory();
    in->FitSlicesY();
    auto mean = (TH1*) current_dir->Get(basename+"_1");
    auto sigma = (TH1*) current_dir->Get(basename+"_2");
    auto out = (TH1 *) mean->Clone(basename+"res");
    for(int ix=1; ix<out->GetNbinsX();++ix) {
        out->SetBinError(ix, sigma->GetBinContent(ix));
    }
    return out;
}

bool match_e(const TLorentzVector &pa, const TLorentzVector &pb, const TLorentzVector &ma, const TLorentzVector &mb)
{
    // pos + direction
    // if (fabs(pa.X() - pb.X() + 0.1) > 0.5) return false;
    // if (fabs(pa.Y() - pb.Y()) > 1.) return false;
    // if (fabs(pa.Z() - pb.Z()) > 1.) return false;
    // if (fabs(ma.Theta() - mb.Theta()) > 0.5) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 0.5) return false;

    // dist
    auto dist = (pa-pb).Vect();
    if (dist.Mag()>1.) return false;

    // if (fabs(ma.Theta() - mb.Theta()) > 1.0) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 1.0) return false;

    return true;
}

void leading_e(
const char *input = "data/nue_overlay_run2.root"
// const char *input = "data/nu_overlay_run2.root"
)
{
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    auto *tf = TFile::Open(input, "read");
    auto *dir = (TDirectoryFile *) tf->Get("wcpselection");
    // auto *T_PFDump = (TTree *) dir->Get("T_PFeval");
    auto *T_PFDump = (TTree *) dir->Get("T_PFDump");
    auto *T_PFeval = (TTree *) dir->Get("T_PFeval");
    T_PFDump->AddFriend(T_PFeval);

    const int target_pdg = 11; // 11, 13, 2212

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

    float truth_corr_nuvtxX;
    float truth_corr_nuvtxY;
    float truth_corr_nuvtxZ;
    float reco_nuvtxX;
    float reco_nuvtxY;
    float reco_nuvtxZ;

    T_PFDump->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_PFDump->SetBranchAddress("truth_pdg", &truth_pdg);
    T_PFDump->SetBranchAddress("truth_mother", &truth_mother);
    T_PFDump->SetBranchAddress("truth_startXYZT", &truth_startXYZT);
    T_PFDump->SetBranchAddress("truth_startMomentum", &truth_startMomentum);
    T_PFDump->SetBranchAddress("truth_daughters", &truth_daughters);

    T_PFDump->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_PFDump->SetBranchAddress("reco_pdg", &reco_pdg);
    T_PFDump->SetBranchAddress("reco_mother", &reco_mother);
    T_PFDump->SetBranchAddress("reco_startXYZT", &reco_startXYZT);
    T_PFDump->SetBranchAddress("reco_startMomentum", &reco_startMomentum);
    T_PFDump->SetBranchAddress("reco_daughters", &reco_daughters);

    T_PFDump->SetBranchAddress("truth_corr_nuvtxX", &truth_corr_nuvtxX);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxY", &truth_corr_nuvtxY);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxZ", &truth_corr_nuvtxZ);
    T_PFDump->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);
    T_PFDump->SetBranchAddress("reco_nuvtxY", &reco_nuvtxY);
    T_PFDump->SetBranchAddress("reco_nuvtxZ", &reco_nuvtxZ);

    TH1F *h_ndaughter_truth = new TH1F("h_ndaughter_truth", "h_ndaughter_truth", 100, -0.5, 99.5);
    TH1F *h_ndaughter_reco = new TH1F("h_ndaughter_reco", "h_ndaughter_reco", 100, -0.5, 99.5);

    TH1F *h_reco_m_truth = new TH1F("h_reco_m_truth", "h_reco_m_truth", 200, -3, 3);
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 4, 100, -1, 1); // energy
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 100, -TMath::Pi(), TMath::Pi()); // theta
    TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 200, -3, 3); // pos

    TH1F *h_truth_e_all = new TH1F("h_truth_e_all", "h_truth_e_all", 100, 0, 3);
    TH1F *h_truth_e_match = new TH1F("h_truth_e_match", "h_truth_e_match", 100, 0, 3);

    TH2F *h_dtheta_theta =
        new TH2F("h_dtheta_theta", "h_dtheta_theta", 200, -TMath::Pi(), TMath::Pi(), 200, -TMath::Pi(), TMath::Pi());

    int counter_has_truth_e = 0;
    int counter_has_reco_e = 0;
    // for (int ientry = 0; ientry < T_PFDump->GetEntries(); ++ientry) {
    for (int ientry = 0; ientry < 100000; ++ientry) {
        T_PFDump->GetEntry(ientry);
        if (ientry % 1000 == 0) cout << "processing: " << ientry / 10000. * 100 << "%" << endl;

        TVector3 truth_nuvtx(truth_corr_nuvtxX, truth_corr_nuvtxY, truth_corr_nuvtxZ);
        TVector3 reco_nuvtx(reco_nuvtxX, reco_nuvtxY, reco_nuvtxZ);
        // truth_nuvtx.Print();
        // reco_nuvtx.Print();
        if((reco_nuvtx-truth_nuvtx).Mag() > 1.0) continue;

        TLorentzVector pos_truth;
        TLorentzVector mom_truth;
        float current_max_energy_truth = 0;
        for (int itruth = 0; itruth < truth_Ntrack; ++itruth) {
            h_ndaughter_truth->Fill(truth_daughters->at(itruth).size());
            if (truth_mother[itruth] != 0) continue;
            if (truth_pdg[itruth] != target_pdg) continue;
            if (truth_startMomentum[itruth][3] < current_max_energy_truth) continue;
            current_max_energy_truth = truth_startMomentum[itruth][3];
            pos_truth.SetXYZT(truth_startXYZT[itruth][0], truth_startXYZT[itruth][1], truth_startXYZT[itruth][2],
                              truth_startXYZT[itruth][3]);
            mom_truth.SetXYZT(truth_startMomentum[itruth][0], truth_startMomentum[itruth][1],
                              truth_startMomentum[itruth][2], truth_startMomentum[itruth][3]);
        }

        if (!fv_shwr_vtx(pos_truth)) continue;

        if (current_max_energy_truth > 0) {
            ++counter_has_truth_e;
        }

        TLorentzVector pos_reco;
        TLorentzVector mom_reco;
        float current_max_energy_reco = 0;
        for (int ireco = 0; ireco < reco_Ntrack; ++ireco) {
            h_ndaughter_reco->Fill(reco_daughters->at(ireco).size());
            if (reco_mother[ireco] != 0) continue;
            if (reco_pdg[ireco] != target_pdg) continue;
            if (reco_startMomentum[ireco][3] < current_max_energy_reco) continue;
            current_max_energy_reco = reco_startMomentum[ireco][3];
            pos_reco.SetXYZT(reco_startXYZT[ireco][0], reco_startXYZT[ireco][1], reco_startXYZT[ireco][2],
                             reco_startXYZT[ireco][3]);
            mom_reco.SetXYZT(reco_startMomentum[ireco][0], reco_startMomentum[ireco][1], reco_startMomentum[ireco][2],
                             reco_startMomentum[ireco][3]);
        }
        if (current_max_energy_reco > 0) {
            ++counter_has_reco_e;
        }

        if (current_max_energy_truth > 0) {
            h_truth_e_all->Fill(mom_truth.E());
        }

        if (current_max_energy_truth == 0 || current_max_energy_reco == 0) continue;

        if (!match_e(pos_reco, pos_truth, mom_reco, mom_truth)) continue;
        h_truth_e_match->Fill(mom_truth.E());
        // h_reco_m_truth->Fill(mom_reco.E() - mom_truth.E());
        h_reco_m_truth->Fill(pos_reco.X() - pos_truth.X());
        // h_reco_v_truth->Fill(mom_truth.E(), mom_reco.E());
        // h_reco_v_truth->Fill(mom_truth.E(),(mom_reco.E()-mom_truth.E())/mom_truth.E()); // energy
        // h_reco_v_truth->Fill(mom_truth.E(),mom_reco.Theta()-mom_truth.Theta()); // angle
        h_reco_v_truth->Fill(pos_truth.Y(), pos_reco.Y()); // pos

        // angle
        auto dtheta =
            TMath::ACos(mom_truth.Vect().Dot(mom_reco.Vect()) / mom_truth.Vect().Mag() / mom_reco.Vect().Mag());
        TLorentzVector mom_truth_yzx;
        dtheta = mom_reco.Phi() - mom_truth.Phi();
        mom_truth_yzx.SetXYZT(mom_truth.Z(), mom_truth.Y(), mom_truth.X(), mom_truth.T());
        h_dtheta_theta->Fill(TMath::Pi() / 2 - mom_truth_yzx.Theta(), TMath::Pi() / 2 - mom_truth.Theta());
    }
    cout << "truth e: " << counter_has_truth_e << "; reco e: " << counter_has_reco_e << endl;

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
    // h_reco_v_truth->SetTitle(";Z^{truth} [cm];Z^{reco}-Z^{truth} [cm]");
    // h_reco_v_truth->SetTitle(";#phi^{truth} [rad];#phi^{reco}-#phi^{truth} [rad]");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco}-E^{truth} [GeV]");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco} [GeV]");
    h_reco_v_truth->SetTitle(";E^{truth} [GeV];#Delta #theta");
    // h_reco_v_truth->SetStats(0);
    h_reco_v_truth->Draw("colz");
    // auto current_dir = h_reco_v_truth->GetDirectory();
    // h_reco_v_truth->FitSlicesY(0, 0, -1, 20);
    // auto h_reco_v_truth_1 = (TH1*) current_dir->Get("h_reco_v_truth_1");
    auto h_reco_v_truth_1 = resolution_hist(h_reco_v_truth);
    h_reco_v_truth_1->Draw("e,same");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->SetGrid();
    h_reco_m_truth->SetTitle(";X^{reco}-X^{truth} [cm]");
    // h_reco_m_truth->SetTitle(";#phi^{reco}-#phi^{truth} [rad]");
    // h_reco_m_truth->SetTitle(";E^{reco}-E^{truth} [GeV]");
    // h_reco_m_truth->SetStats(0);
    h_reco_m_truth->Draw();

    // TCanvas *c4 = new TCanvas("c4", "c4");
    // h_ndaughter_truth->Draw();
    // TCanvas *c5 = new TCanvas("c5", "c5");
    // h_ndaughter_reco->Draw();

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->SetGrid();
    h_dtheta_theta->SetTitle("; #theta_{YZ}; #theta_{XY}");
    // h_reco_m_truth->SetStats(0);
    h_dtheta_theta->Draw("colz");
}
