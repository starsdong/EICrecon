#include <iostream>
// #include <fstream>

// #include <string>
#include <vector>
// #include <sstream>
// #include <cstdlib>
using namespace std;

void plot_hists(){

    //Define Style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);

    //Get ROOT file
    TFile *f = new TFile("eicrecon.root");
    TH2 *h1a = (TH2*) f->Get("track_qa/h1a");
    TH1 *hchi2 = (TH1*) f->Get("track_qa/hchi2");
    TH1 *heta = (TH1*) f->Get("track_qa/heta");
    TH1 *hp = (TH1*) f->Get("track_qa/hp");
    TH1 *hpt = (TH1*) f->Get("track_qa/hpt");
    TH1 *hhits = (TH1*) f->Get("track_qa/hhits");
    TH1 *hNDF = (TH1*) f->Get("track_qa/hNDF");
    TH1 *hchi2_by_hits = (TH1*) f->Get("track_qa/hchi2_by_hits");
    TH1 *hchi2_by_NDF = (TH1*) f->Get("track_qa/hchi2_by_NDF");

    TH2 *hchi2_vs_eta = (TH2*) f->Get("track_qa/hchi2_vs_eta");
    TH2 *hchi2_vs_hits = (TH2*) f->Get("track_qa/hchi2_vs_hits");
    TH2 *hchi2_vs_hits_zoomed = (TH2*) f->Get("track_qa/hchi2_vs_hits_zoomed");
    // vector<TH2*> hchi2_vs_hits_etabins = (vector<TH2*>) f->Get("track_qa/hchi2_vs_hits_etabins");
    TH2 *hhits_vs_eta = (TH2*) f->Get("track_qa/hhits_vs_eta");
    TH2 *htracks_vs_eta = (TH2*) f->Get("track_qa/htracks_vs_eta");
    TH3 *heta_vs_p_vs_chi2 = (TH3*) f->Get("track_qa/heta_vs_p_vs_chi2");
    TH2 *hmeasptrack_vs_eta = (TH2*) f->Get("track_qa/hmeasptrack_vs_eta");
    TH2 *hmeasptrack_vs_hits = (TH2*) f->Get("track_qa/hmeasptrack_vs_hits");
    TH2 *hmeasptrack_vs_chi2perNDF = (TH2*) f->Get("track_qa/hmeasptrack_vs_chi2perNDF");
    TH2 *hmeasptrack_vs_calstates = (TH2*) f->Get("track_qa/hmeasptrack_vs_calstates");

    TH2 *hmeaschi2_vs_chi2 = (TH2*) f->Get("track_qa/hmeaschi2_vs_chi2");
    TH2 *hmeaschi2_vs_eta = (TH2*) f->Get("track_qa/hmeaschi2_vs_eta");
    TH2 *hmeaschi2_vs_hits = (TH2*) f->Get("track_qa/hmeaschi2_vs_hits");

       //use 0.5i-4 to get lowerbound and 0.5i-3.5 to get upper bound
    vector<TH2*> hchi2_vs_hits_etabins(16);
    vector<TH2*> hmeasptrack_vs_hits_etabins(16);
    vector<TH2*> hmeasptrack_vs_hits_etabins_zoomed(16);
    vector<TH2*> hmeasptrack_vs_chi2perNDF_etabins(16);
    for (int i=0; i<16; i++){
        hchi2_vs_hits_etabins[i] = (TH2*) f->Get(
            TString::Format("track_qa/eta_bins/hchi2_vs_hits_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5));
        hmeasptrack_vs_hits_etabins[i] = (TH2*) f->Get(
            TString::Format("track_qa/eta_bins/hmeasptrack_vs_hits_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5));
        hmeasptrack_vs_hits_etabins_zoomed[i] = (TH2*) f->Get(
            TString::Format("track_qa/eta_bins/hmeasptrack_vs_hits_eta_zoomed_%.1f_%.1f", 0.5*i-4,0.5*i-3.5));
        hmeasptrack_vs_chi2perNDF_etabins[i] = (TH2*) f->Get(
            TString::Format("track_qa/eta_bins/hmeasptrack_vs_chi2perNDF_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5));
    }

    TH2 *hholes_vs_hits = (TH2*) f->Get("track_qa/hholes_vs_hits");
    TH2 *houtliers_vs_hits = (TH2*) f->Get("track_qa/houtliers_vs_hits");
    TH2 *hsummation = (TH2*) f->Get("track_qa/hsummation");
    TH2 *hsummation2 = (TH2*) f->Get("track_qa/hsummation2");

    

    //Make plots
    TCanvas *c1a = new TCanvas("c1a");
    h1a->Draw("colz");

    TPaveText* tex_gen = new TPaveText(0.2,0.7,0.5,0.85,"NDCNB");
    tex_gen->AddText("Single Electrons generated:");
    tex_gen->AddText("1 GeV < E < 10 GeV");
    // tex_gen->AddText("#theta = 90^{o}, 0^{o} < #phi < 360^{o}");
    tex_gen->AddText("-4 < #eta < 4, 0^{o} < #phi < 360^{o}");
    tex_gen->SetFillStyle(4000);tex_gen->SetTextFont(63);tex_gen->SetTextSize(20);
    tex_gen->Draw();

    TPaveText* tex_zoom = new TPaveText(0.7,0.9,0.8,0.95,"NDCNB");
    tex_zoom->AddText("Zoomed in");
    tex_zoom->SetFillStyle(4000);tex_zoom->SetTextFont(63);tex_zoom->SetTextSize(20);tex_zoom->SetTextColor(kRed);

    TCanvas *c2a = new TCanvas("c2a");
    hchi2->Draw();
    tex_gen->Draw();

    TCanvas *c2b = new TCanvas("c2b");
    c2b->Divide(2,2);
    c2b->cd(1);
    hhits->Draw();
    c2b->cd(2);
    hNDF->Draw();
    c2b->cd(3);
    hchi2_by_hits->Draw();
    c2b->cd(4);
    hchi2_by_NDF->Draw();

    TCanvas *c2c = new TCanvas("c2c","c2c",1200,900);
    c2c->Divide(2,2);
    c2c->cd(1);
    heta->Draw();
    c2c->cd(2);
    hp->Draw();
    c2c->cd(3);
    hpt->Draw();

    TCanvas *c3a = new TCanvas("c3a");
    hchi2_vs_eta->Draw("colz");


    vector<TPaveText*> tex_eta(16);
    for (int i=0; i<16; i++){
        tex_eta[i] = new TPaveText(0.5,0.7,0.7,0.85,"NDCNB");
        tex_eta[i]->AddText(TString::Format("%.1f < #eta < %.1f", 0.5*i-4, 0.5*i-3.5));
        tex_eta[i]->SetFillStyle(4000);tex_eta[i]->SetTextFont(63);tex_eta[i]->SetTextSize(10);
    }

    TCanvas *c4a = new TCanvas("c4a");
    hchi2_vs_hits->Draw("colz");
    TCanvas *c4b = new TCanvas("c4b");
    hchi2_vs_hits_zoomed->Draw("colz");
    TCanvas *c4c = new TCanvas("c4c");
    c4c->Divide(4,4);
    for (int i=0; i<16; i++){
        c4c->cd(i+1);
        hchi2_vs_hits_etabins[i]->Draw("colz");
        tex_eta[i]->Draw();
    }

    TCanvas *c5a = new TCanvas("c5a");
    hhits_vs_eta->Draw("colz");

    TCanvas *c6a = new TCanvas("c6a");
    htracks_vs_eta->Draw("colz");

    TCanvas *c7a = new TCanvas("c7a");
    heta_vs_p_vs_chi2->Draw("lego2");

    TCanvas *c8a = new TCanvas("c8a");
    hmeasptrack_vs_eta->Draw("colz");


    auto fdiagline = new TF1("fdiagline","x",0,10);
    fdiagline->SetLineWidth(3);fdiagline->SetLineColor(kRed);
    auto fdiagline_thin = new TF1("fdiagline_thin","x",0,10);
    fdiagline_thin->SetLineWidth(2);fdiagline_thin->SetLineColor(kRed);

    TCanvas *c9a = new TCanvas("c9a");
    hmeasptrack_vs_hits->Draw("colz");
    fdiagline->Draw("same");
    TCanvas *c9b = new TCanvas("c9b");
    c9b->Divide(4,4);
    for (int i=0; i<16; i++){
        c9b->cd(i+1);
        hmeasptrack_vs_hits_etabins[i]->Draw("colz");
        tex_eta[i]->Draw();
        fdiagline_thin->Draw("same");
    }
    TCanvas *c9c = new TCanvas("c9c");
    c9c->Divide(4,4);
    for (int i=0; i<16; i++){
        c9c->cd(i+1);
        hmeasptrack_vs_hits_etabins_zoomed[i]->Draw("colz");
        tex_eta[i]->Draw();
        fdiagline_thin->Draw("same");
    }

    TCanvas *c10a = new TCanvas("c10a");
    hmeasptrack_vs_chi2perNDF->Draw("colz");
    TCanvas *c10b = new TCanvas("c10b");
    c10b->Divide(4,4);
    for (int i=0; i<16; i++){
        c10b->cd(i+1);
        hmeasptrack_vs_chi2perNDF_etabins[i]->Draw("colz");
        tex_eta[i]->Draw();
    }
    
    TCanvas *c11a = new TCanvas("c11a");
    hmeasptrack_vs_calstates->Draw("colz");



    TCanvas *c12a = new TCanvas("c12a");
    hmeaschi2_vs_chi2->Draw("colz");

    TCanvas *c13a = new TCanvas("c13a");
    hmeaschi2_vs_eta->Draw("colz");

    TCanvas *c14a = new TCanvas("c14a");
    hmeaschi2_vs_hits->Draw("colz");

    TCanvas *c15a = new TCanvas("c15a");
    hholes_vs_hits->Draw("colz");

    TCanvas *c16a = new TCanvas("c16a");
    houtliers_vs_hits->Draw("colz");

    TCanvas *c17a = new TCanvas("c17a");
    hsummation->Draw("colz");

    TCanvas *c18a = new TCanvas("c18a");
    hsummation2->Draw("colz");


    //Print plots to file
    c1a->Print("plot_hists_etarange_flat.pdf[");
    c1a->Print("plot_hists_etarange_flat.pdf");
    c2a->Print("plot_hists_etarange_flat.pdf");
    c2b->Print("plot_hists_etarange_flat.pdf");
    c2c->Print("plot_hists_etarange_flat.pdf");
    c3a->Print("plot_hists_etarange_flat.pdf");
    c4a->Print("plot_hists_etarange_flat.pdf");
    c4b->Print("plot_hists_etarange_flat.pdf");
    c4c->Print("plot_hists_etarange_flat.pdf");
    c5a->Print("plot_hists_etarange_flat.pdf");
    c6a->Print("plot_hists_etarange_flat.pdf");
    c7a->Print("plot_hists_etarange_flat.pdf");
    c8a->Print("plot_hists_etarange_flat.pdf");
    c9a->Print("plot_hists_etarange_flat.pdf");
    c9b->Print("plot_hists_etarange_flat.pdf");
    c9c->Print("plot_hists_etarange_flat.pdf");
    c10a->Print("plot_hists_etarange_flat.pdf");
    c10b->Print("plot_hists_etarange_flat.pdf");
    c11a->Print("plot_hists_etarange_flat.pdf");
    c12a->Print("plot_hists_etarange_flat.pdf");
    c13a->Print("plot_hists_etarange_flat.pdf");
    c14a->Print("plot_hists_etarange_flat.pdf");
    c15a->Print("plot_hists_etarange_flat.pdf");
    c16a->Print("plot_hists_etarange_flat.pdf");
    c17a->Print("plot_hists_etarange_flat.pdf");
    c18a->Print("plot_hists_etarange_flat.pdf");
    c18a->Print("plot_hists_etarange_flat.pdf]");
}