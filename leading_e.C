#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>

using namespace std;

enum LIMITS {
    MAX_TRACKS = 30000,
};

bool fv_shwr_vtx (const TLorentzVector & pos) {
    if(pos.X() > 23 && pos.X() < 233 &&
            pos.Y() > -93 && pos.Y() < 93 &&
            pos.Z() > 23 && pos.Z() < 1014
      ) {
        return true;
    }

    return false;
}

bool match_e(
        const TLorentzVector & pa, const TLorentzVector &pb,
        const TLorentzVector & ma, const TLorentzVector &mb
        ) {
    // pos + direction
    if (fabs(pa.X()-pb.X()+0.1)>0.5) return false;
    if (fabs(pa.Y()-pb.Y())>1.) return false;
    if (fabs(pa.Z()-pb.Z())>1.) return false;
    if (fabs(ma.Theta()-mb.Theta())>0.5) return false;
    if (fabs(ma.Phi()-mb.Phi())>0.5) return false;

    // dist
    //auto dist = (pa-pb).Vect();
    //if (dist.Mag()>1.) return false;

    return true;
}

void leading_e(
        const char* input = "nue_overlay_run2.root"
        //const char* input = "nue_overlay_run1_test_3000.root"
        ) {
    auto *tf = TFile::Open(input, "read");
    auto *dir = (TDirectoryFile*) tf->Get("wcpselection");
    auto *T_PFDump = (TTree*) dir->Get("T_PFDump");

    int truth_Ntrack;
    int truth_id[MAX_TRACKS];
    int truth_pdg[MAX_TRACKS];
    int truth_process[MAX_TRACKS];
    int truth_mother[MAX_TRACKS];
    float truth_startXYZT[MAX_TRACKS][4];
    float truth_endXYZT[MAX_TRACKS][4];
    float truth_startMomentum[MAX_TRACKS][4];
    float truth_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *truth_daughters;

    int reco_Ntrack;
    int reco_id[MAX_TRACKS];
    int reco_pdg[MAX_TRACKS];
    int reco_process[MAX_TRACKS];
    int reco_mother[MAX_TRACKS];
    float reco_startXYZT[MAX_TRACKS][4];
    float reco_endXYZT[MAX_TRACKS][4];
    float reco_startMomentum[MAX_TRACKS][4];
    float reco_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *reco_daughters;

    T_PFDump->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_PFDump->SetBranchAddress("truth_pdg", &truth_pdg);
    T_PFDump->SetBranchAddress("truth_mother", &truth_mother);
    T_PFDump->SetBranchAddress("truth_startXYZT", &truth_startXYZT);
    T_PFDump->SetBranchAddress("truth_startMomentum", &truth_startMomentum);

    T_PFDump->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_PFDump->SetBranchAddress("reco_pdg", &reco_pdg);
    T_PFDump->SetBranchAddress("reco_mother", &reco_mother);
    T_PFDump->SetBranchAddress("reco_startXYZT", &reco_startXYZT);
    T_PFDump->SetBranchAddress("reco_startMomentum", &reco_startMomentum);

    TH1F* h_reco_m_truth = new TH1F("h_reco_m_truth","h_reco_m_truth",200,-3,3);
    TH2F* h_reco_v_truth = new TH2F("h_reco_v_truth","h_reco_v_truth",200,0,4,200,0,4);

    TH1F* h_truth_e_all = new TH1F("h_truth_e_all","h_truth_e_all",200,0,4);
    TH1F* h_truth_e_match = new TH1F("h_truth_e_match","h_truth_e_match",200,0,4);

    int counter_has_truth_e = 0;
    int counter_has_reco_e = 0;
    for (int ientry=0; ientry<T_PFDump->GetEntries(); ++ientry) {
        //for (int ientry=0; ientry<10000; ++ientry) {
        T_PFDump->GetEntry(ientry);
        if(ientry%1000==0) cout << "processing: " << ientry/10000.*100 << "%" << endl;

        TLorentzVector pos_truth;
        TLorentzVector mom_truth;
        float current_max_energy_truth = 0;
        for(int itruth=0; itruth<truth_Ntrack; ++itruth) {
            if(truth_mother[itruth]!=0) continue;
            if(truth_pdg[itruth]!=11) continue;
            if(truth_startMomentum[itruth][3]<current_max_energy_truth) continue;
            current_max_energy_truth = truth_startMomentum[itruth][3];
            pos_truth.SetXYZT(
                    truth_startXYZT[itruth][0],
                    truth_startXYZT[itruth][1],
                    truth_startXYZT[itruth][2],
                    truth_startXYZT[itruth][3]
                    );
            mom_truth.SetXYZT(
                    truth_startMomentum[itruth][0],
                    truth_startMomentum[itruth][1],
                    truth_startMomentum[itruth][2],
                    truth_startMomentum[itruth][3]
                    );
        }

        if(!fv_shwr_vtx(pos_truth)) continue;

        if(current_max_energy_truth>0) {
            ++counter_has_truth_e;
        }

        TLorentzVector pos_reco;
        TLorentzVector mom_reco;
        float current_max_energy_reco = 0;
        for(int ireco=0; ireco<reco_Ntrack; ++ireco) {
            if(reco_mother[ireco]!=0) continue;
            if(reco_pdg[ireco]!=11) continue;
            if(reco_startMomentum[ireco][3]<current_max_energy_reco) continue;
            current_max_energy_reco = reco_startMomentum[ireco][3];
            pos_reco.SetXYZT(
                    reco_startXYZT[ireco][0],
                    reco_startXYZT[ireco][1],
                    reco_startXYZT[ireco][2],
                    reco_startXYZT[ireco][3]
                    );
            mom_reco.SetXYZT(
                    reco_startMomentum[ireco][0],
                    reco_startMomentum[ireco][1],
                    reco_startMomentum[ireco][2],
                    reco_startMomentum[ireco][3]
                    );
        }
        if(current_max_energy_reco>0) {
            ++counter_has_reco_e;
        }

        if(current_max_energy_truth>0) {
            h_truth_e_all->Fill(mom_truth.E());
        }

        if(current_max_energy_truth==0 || current_max_energy_reco==0) continue;

        if(!match_e(pos_reco,pos_truth,mom_reco,mom_truth)) continue;
        h_truth_e_match->Fill(mom_truth.E());
        h_reco_m_truth->Fill(mom_reco.E()-mom_truth.E());
        h_reco_v_truth->Fill(mom_truth.E(),mom_reco.E());
        //h_reco_v_truth->Fill(mom_truth.E(),mom_reco.E()-mom_truth.E());
    }
    cout << "truth e: " << counter_has_truth_e << "; reco e: " << counter_has_reco_e << endl;

    const int LINE_WIDTH = 2;

    TCanvas *c0 = new TCanvas("c0","c0");
    h_truth_e_all->SetLineColor(kBlack);
    h_truth_e_all->SetLineWidth(LINE_WIDTH);
    h_truth_e_all->SetTitle(";E^{truth} [GeV]");
    h_truth_e_all->Draw();
    h_truth_e_match->SetLineColor(kRed);
    h_truth_e_match->SetLineWidth(LINE_WIDTH);
    h_truth_e_match->Draw("same");

    auto h_truth_e_ratio = (TH1F*)h_truth_e_match->Clone("h_truth_e_ratio");
    h_truth_e_ratio->Divide(h_truth_e_all);
    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetGrid();
    h_truth_e_ratio->SetLineColor(kRed);
    h_truth_e_ratio->SetLineWidth(LINE_WIDTH);
    h_truth_e_ratio->SetTitle(";E^{truth} [GeV];Efficiency");
    h_truth_e_ratio->SetStats(0);
    h_truth_e_ratio->Draw("");


    TCanvas *c2 = new TCanvas("c2","c2");
    c2->SetGrid();
    //h_reco_v_truth->SetTitle(";Z^{truth} [cm];Z^{reco}-Z^{truth} [cm]");
    //h_reco_v_truth->SetTitle(";#phi^{truth} [rad];#phi^{reco}-#phi^{truth} [rad]");
    //h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco}-E^{truth} [GeV]");
    h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco} [GeV]");
    //h_reco_v_truth->SetStats(0);
    h_reco_v_truth->Draw("colz");

    TCanvas *c3 = new TCanvas("c3","c3");
    c3->SetGrid();
    //h_reco_m_truth->SetTitle(";Z^{reco}-Z^{truth} [cm]");
    //h_reco_m_truth->SetTitle(";#phi^{reco}-#phi^{truth} [rad]");
    h_reco_m_truth->SetTitle(";E^{reco}-E^{truth} [GeV]");
    //h_reco_m_truth->SetStats(0);
    h_reco_m_truth->Draw();
    }
