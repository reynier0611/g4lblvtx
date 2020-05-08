/*
Winston DeGraw (wdegraw@lbl.gov)
Rey Cruz-Torres (reynier@lbl.gov)
*/

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker );
// ============================================================================================================================================
void analysis_resolution(){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile("muons_AllSiBarrel_FastTrackingEval.root");
	TTree * T = (TTree*) F -> Get("tracks");
	float gpx, gpy, gpz, px, py, pz;
	T -> SetBranchAddress("gpx",&gpx);
	T -> SetBranchAddress("gpy",&gpy);
	T -> SetBranchAddress("gpz",&gpz);
	T -> SetBranchAddress("px" ,&px );
	T -> SetBranchAddress("py" ,&py );
	T -> SetBranchAddress("pz" ,&pz );
	int nEntries = T -> GetEntries();
	// -------------------------------------------------------------
	// Defining binning that will subsequently be used
	float eta_binning[] = {0.,0.2,0.4,0.6,0.8,1.1};
	const int size_eta_binning = sizeof(eta_binning)/sizeof(*eta_binning);
	float p_binning[] = {0.,0.5,1.,1.5,2.,3.,5.,7.,9.,12.,15.,20.,25.};
	const int size_p_binning = sizeof(p_binning)/sizeof(*p_binning);
	// -------------------------------------------------------------
	// Defining histograms
	TH1F *** h1_dp_p = new TH1F**[size_eta_binning-1];
	TH1F ** h1_resolution = new TH1F*[size_eta_binning-1];
	for(int et = 0 ; et < size_eta_binning-1 ; et++){
		h1_dp_p[et] = new TH1F*[size_p_binning-1];
		for(int p = 0 ; p < size_p_binning-1 ; p++){
			h1_dp_p[et][p] = new TH1F(Form("h1_dp_p_%i_%i",et,p),Form("%.1f < #eta < %.1f, %.1f < p < %.1f GeV/c;#Deltap/p;Counts",eta_binning[et],eta_binning[et+1],p_binning[p],p_binning[p+1]),100,-.15,.15);
		}

		h1_resolution[et] = new TH1F(Form("h1_resolution_%i",et),";p [GeV/c];#Deltap/p",size_p_binning-1,p_binning);
	}
	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width[size_eta_binning-1][size_p_binning-1] = {{0}};
	float error[size_eta_binning-1][size_p_binning-1] = {{0}};
	TString label[size_eta_binning-1];
	TF1 *** f_gaus = new TF1**[size_eta_binning-1];
        for(int et = 0 ; et < size_eta_binning-1 ; et++){
                f_gaus[et] = new TF1*[size_p_binning-1];
                for(int p = 0 ; p < size_p_binning-1 ; p++){
                        f_gaus[et][p] = new TF1(Form("f_gaus_%i_%i",et,p),"gaus",-.05,.05);
                }
		label[et] = Form("%.1f < |#eta| < %.1f",eta_binning[et],eta_binning[et+1]);
        }
	// -------------------------------------------------------------
	// Loop over entries of the tree
	for(int ev = 0 ; ev < nEntries ; ev++){
		T -> GetEntry(ev);
		if(ev%10000==0) cout << "Looping over entry " << ev << " out of " << nEntries << endl;

		float gtheta = TMath::ACos(gpz/sqrt(gpx*gpx+gpy*gpy+gpz*gpz));
		float geta = -TMath::Log(TMath::Tan(gtheta/2.));
		float p_reco = sqrt(px*px+py*py+pz*pz);
		float p_truth = sqrt(gpx*gpx+gpy*gpy+gpz*gpz);
		float dp_p = (p_reco-p_truth)/p_truth;

		for(int et = 0 ; et < size_eta_binning-1 ; et++){
			if( abs(geta) >  eta_binning[et] &&  abs(geta) <= eta_binning[et+1] ){
			
				for(int p = 0 ; p < size_p_binning-1 ; p++){
					if( p_truth > p_binning[p] && p_truth <= p_binning[p+1] ){
						h1_dp_p[et][p] -> Fill(dp_p);
					}	
				}
			}
		}
	}
	// -------------------------------------------------------------
	// Plotting histograms
	TCanvas ** c_fits = new TCanvas*[size_eta_binning-1];

	int color[] = {1,2,62,8,93};
	int marker[] = {20,21,23,24,25};

	for(int et = 0 ; et < size_eta_binning-1 ; et++){
		c_fits[et] = new TCanvas(Form("c_fits_%i",et),Form("%.1f<eta<%.1f",eta_binning[et],eta_binning[et+1]),1000,800);
		c_fits[et] -> Divide(4,3);
		for(int p = 0 ; p < size_p_binning-1 ; p++){
			c_fits[et] -> cd(p+1);
			h1_dp_p[et][p] -> Draw();								// Drawing previously filled histogram
			h1_dp_p[et][p] -> Fit(Form("f_gaus_%i_%i",et,p),"R");					// Fitting a gaussian to the histogram
			width[et][p] = f_gaus[et][p] -> GetParameter(2);					// Saving the standard deviation in a variable
			float fit_err = f_gaus[et][p] -> GetParError(2);			// Saving standard deviation fit error in a variable
			float chi2_dof = (f_gaus[et][p] -> GetChisquare())/(f_gaus[et][p] -> GetNDF()); // Saving chi^2 per dof in a variable
			error[et][p] = fit_err*chi2_dof;		// Defining error on width as standard deviation x chi^2/dof

			h1_resolution[et] -> SetBinContent(p+1,width[et][p]);					// Saving resolution value in dp/p vs. p histogram
			h1_resolution[et] -> SetBinError(p+1,error[et][p]);					// Saving error
		}
		prettyTH1F( h1_resolution[et] , color[et] , marker[et] );					// Editing the histogram
	}

	TCanvas * c1 = new TCanvas("c1");
	h1_resolution[1] -> Draw();
	for(int et = 0 ; et < size_eta_binning-1 ; et++){
		h1_resolution[et] -> Draw("same");
	}

	TLegend * leg1 = new TLegend(0.20,0.6,0.40,0.85);
	leg1 -> SetLineColor(0);
	for(int et = 0 ; et < size_eta_binning-1 ; et++) leg1 -> AddEntry( h1_resolution[et],label[et] );
	leg1 -> Draw("same");
    c1 -> Print("muons_AllSiBarrel_dppVp.pdf");
    
    
    TString out_pdf_name = "muons_AllSiBarrel_fits.pdf";
    for(int et = 0 ; et < size_eta_binning-1 ; et++){
        TString fname = out_pdf_name;
        if(et == 0) fname+="(";
        else if(et == size_eta_binning-2) fname+=")";
        c_fits[et] -> Print(fname);
    }


}


// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker ){
	h1 -> SetLineWidth(2);
	h1 -> SetLineColor(color);
	h1 -> SetMarkerStyle(marker);
	h1 -> SetMarkerColor(color);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
    h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks
}
