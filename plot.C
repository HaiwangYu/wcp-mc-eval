#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>

using namespace std;

enum LIMITS {
    MAX_TRACKS = 30000,
};

bool tlv_match(const TLorentzVector & a, const TLorentzVector & b) {
    if (fabs(a.X()-b.X())>0.5) return false;
    if (fabs(a.Y()-b.Y())>0.5) return false;
    if (fabs(a.Z()-b.Z())>0.5) return false;
    if (fabs(a.T()-b.T())>1.0) return false;
    return true;
}

void plot(
        const char* input = "reco2_checkout_hist.root"
        ) {
    auto *tf = TFile::Open(input, "read");
    auto *dir = (TDirectoryFile*) tf->Get("wcpselection");
    auto *T_MCEval = (TTree*) dir->Get("T_MCEval");

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

    T_MCEval->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_MCEval->SetBranchAddress("truth_pdg", &truth_pdg);
    T_MCEval->SetBranchAddress("truth_mother", &truth_mother);
    T_MCEval->SetBranchAddress("truth_startMomentum", &truth_startMomentum);

    T_MCEval->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_MCEval->SetBranchAddress("reco_pdg", &reco_pdg);
    T_MCEval->SetBranchAddress("reco_mother", &reco_mother);
    T_MCEval->SetBranchAddress("reco_startMomentum", &reco_startMomentum);

    const float cut_energy = 0.1;
    const float plot_max_energy = 3.0;

    TH1F* h_truth_energy = new TH1F("h_truth_energy","h_truth_energy",100,0,plot_max_energy);
    TH1F* h_reco_energy = new TH1F("h_reco_energy","h_reco_energy",100,0,plot_max_energy);
    TH1F* h_reco_energy_match = new TH1F("h_reco_energy_match","h_reco_energy_match",100,0,plot_max_energy);

    TH2F* h_truth_energy_reco_energy = new TH2F("h_truth_energy_reco_energy","h_truth_energy_reco_energy",100,0,plot_max_energy,100,0,plot_max_energy);
    TH2F* h_truth_energy_reco_energy_match = new TH2F("h_truth_energy_reco_energy_match","h_truth_energy_reco_energy_match",100,0,plot_max_energy,100,0,plot_max_energy);

    for (int ientry=0; ientry<T_MCEval->GetEntries(); ++ientry) {
    //for (int ientry=0; ientry<100; ++ientry) {
        T_MCEval->GetEntry(ientry);
        cout << "ientry: " << ientry << endl;
        for(int itruth=0; itruth<truth_Ntrack; ++itruth) {
            if(truth_mother[itruth]!=0) continue;
            if(truth_pdg[itruth]!=11) continue;
            TLorentzVector tlv_truth(
                    truth_startMomentum[itruth][0],
                    truth_startMomentum[itruth][1],
                    truth_startMomentum[itruth][2],
                    truth_startMomentum[itruth][3]
                    );
            h_truth_energy->Fill(tlv_truth.E());

            for(int ireco=0; ireco<truth_Ntrack; ++ireco) {
                if(reco_mother[ireco]!=0) continue;
                if(reco_pdg[ireco]!=11) continue;
                TLorentzVector tlv_reco(
                        reco_startMomentum[ireco][0],
                        reco_startMomentum[ireco][1],
                        reco_startMomentum[ireco][2],
                        reco_startMomentum[ireco][3]
                        );
                h_truth_energy_reco_energy->Fill(tlv_truth.E(),tlv_reco.E());
                if(tlv_reco.E()<cut_energy) continue;
                if(!tlv_match(tlv_truth,tlv_reco)) continue;
                h_truth_energy_reco_energy_match->Fill(tlv_truth.E(),tlv_reco.E());
                h_reco_energy_match->Fill(tlv_reco.E());
            }
        }

        for(int ireco=0; ireco<truth_Ntrack; ++ireco) {
            if(reco_mother[ireco]!=0) continue;
            if(reco_pdg[ireco]!=11) continue;
            TLorentzVector tlv_reco(
                    reco_startMomentum[ireco][0],
                    reco_startMomentum[ireco][1],
                    reco_startMomentum[ireco][2],
                    reco_startMomentum[ireco][3]
                    );
            if(tlv_reco.E()<cut_energy) continue;
            h_reco_energy->Fill(tlv_reco.E());
        }

    }

    const int LINE_WIDTH = 2;

    TCanvas *c0 = new TCanvas("c0","c0");
    h_reco_energy->SetLineColor(kBlue);
    h_reco_energy->SetLineWidth(LINE_WIDTH);
    h_reco_energy->Draw();
    h_truth_energy->SetLineColor(kBlack);
    h_truth_energy->SetLineWidth(LINE_WIDTH);
    h_truth_energy->Draw("same");
    h_reco_energy_match->SetLineColor(kRed);
    h_reco_energy_match->SetLineWidth(LINE_WIDTH);
    h_reco_energy_match->Draw("same");

    TCanvas *c1 = new TCanvas("c1","c1");
    h_truth_energy_reco_energy->SetTitle(";E_{e}^{truth} [GeV];E_{e}^{reco} [GeV]");
    h_truth_energy_reco_energy->SetStats(0);
    h_truth_energy_reco_energy->Draw("colz");

    TCanvas *c2 = new TCanvas("c2","c2");
    h_truth_energy_reco_energy_match->SetTitle(";E_{e}^{truth} [GeV];E_{e}^{reco} [GeV]");
    h_truth_energy_reco_energy_match->SetStats(0);
    h_truth_energy_reco_energy_match->Draw("colz");
}
