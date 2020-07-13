vector<float> thetavec[5]; vector<float> phivec[5];   vector<float> x0vec_max[5];       vector<float> x0vec_min[5];    vector<float> lamda0vec[5];
int loaded = 0;
float ** thetaarr; float ** etaarr; float ** phiarr; float ** x0arr;    float ** x0arr_copy; float ** x0Earr; float ** lamda0arr;

TGraphErrors ** gr_th = new TGraphErrors*[5];
TGraphErrors ** gr_et = new TGraphErrors*[5];

int color[] = {1,6,2,62,8,13};

TString part[] = {"total","beampipe","barrel","forward","backward"};
const int size_part = sizeof(part)/sizeof(*part);


// Forward-delaring functions
void LoadData(const char *datfile);
void PrepareVariables();
void CreatePlots();
void PrettyTGraphErrors( TGraphErrors * gP , int color , int marker , float xmin , float xmax , float ymin , float ymax , TString xtit , TString ytit , TString tit );
// ========================================================================================================================================
void material_scan(){

	for(int i = 0 ; i < size_part ; i++) LoadData("data/matscan_"+part[i]+"_phiall.dat");
	PrepareVariables();
	CreatePlots();

	// -------------------------------------------------
	TLegend * leg = new TLegend(0.43,0.5,0.62,0.8);
	leg -> SetLineColor(0);
	for(int i = 0 ; i < size_part ; i++) leg -> AddEntry(gr_et[i],part[i]);

	// -------------------------------------------------
	// Plotting results
	TCanvas *c1 = new TCanvas("c1","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_th[0]->SetMaximum(0.3);
	gr_th[0]->Draw("ALE3");
	for(int l = 0 ; l < loaded ; l++) gr_th[l] -> Draw("sameLE3");
	leg -> Draw("same");

	TCanvas *c2 = new TCanvas("c2","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[0]->Draw("ALE3");
	for(int l = 0 ; l < loaded ; l++) gr_et[l] -> Draw("sameLE3");
	leg -> Draw("same");

	// -------------------------------------------------

	TCanvas *c3 = new TCanvas("c3","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[0]->Draw("ALE3");
	leg -> Draw("same");

	TCanvas *c4 = new TCanvas("c4","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[0]->Draw("ALE3");
	gr_et[1]->Draw("SAMELE3");
	leg -> Draw("same");

	TCanvas *c5 = new TCanvas("c5","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[0]->Draw("ALE3");
	gr_et[1]->Draw("SAMELE3");
	gr_et[2]->Draw("SAMELE3");
	leg -> Draw("same");

	TCanvas *c6 = new TCanvas("c6","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[2]->Draw("ALE3");

	TCanvas *c7 = new TCanvas("c7","material scan"); gPad -> SetTopMargin(0.03); gPad -> SetRightMargin(0.03);
	gr_et[3]->Draw("ALE3");

	// -------------------------------------------------
	c1 -> Print("results1.pdf(");
	c3 -> Print("results1.pdf" );
	c4 -> Print("results1.pdf" );
	c5 -> Print("results1.pdf" );
	c2 -> Print("results1.pdf)");

	return;
}
// ========================================================================================================================================
void LoadData(const char *datfile){
	FILE *f = fopen(datfile,"r");
	if (!f){cout << "file " << datfile << " cannot be opened" << endl;      return;}
	float theta;    float phi;      float path;     float x0;       float lamda0;
	float thetamin = 10000 ;        float thetamax = -10000 ;       float phimin    = 10000 ;       float phimax    = -10000 ;
	float x0min    = 100000;        float x0max    = -100000;       float lamda0min = 100000;       float lamda0max = -100000;

	float previous_theta = -999.;
	float previous_x0_max;
	float previous_x0_min;

	int loaded_points = 0;

	while(fscanf(f,"%f %f %f %f %f",&theta,&phi,&path,&x0,&lamda0) != EOF){
		if (theta > thetamax)   thetamax = theta;
		if (theta < thetamin)   thetamin = theta;
		if (phi > phimax)       phimax = phi;
		if (phi < phimin)       phimin = phi;
		if (x0 > x0max)         x0max = x0;
		if (x0 < x0min)         x0min = x0;
		if (lamda0 > lamda0max) lamda0max = lamda0;
		if (lamda0 < lamda0min) lamda0min = lamda0;

		/*
		   if( loaded_points==0 ){ 
		   thetavec[loaded].push_back(theta);
		   phivec[loaded].push_back(phi);
		   x0vec_max[loaded].push_back(x0);
		   lamda0vec[loaded].push_back(lamda0);

		   previous_theta = theta;
		   previous_x0_max = x0;
		   previous_x0_min = x0;
		   }
		   else{
		   */
		if(theta != previous_theta){
			thetavec[loaded].push_back(theta);
			phivec[loaded].push_back(phi);
			x0vec_max[loaded].push_back(x0);
			x0vec_min[loaded].push_back(x0);
			lamda0vec[loaded].push_back(lamda0);

			previous_theta = theta;
			previous_x0_max = x0;
			previous_x0_min = x0;
		}
		else{
			if(x0 > previous_x0_max){
				x0vec_max[loaded][x0vec_max[loaded].size()-1] = x0;
				previous_x0_max = x0;
			}
			if(x0 < previous_x0_min){
				x0vec_min[loaded][x0vec_min[loaded].size()-1] = x0;
				previous_x0_min = x0;
			}
		}
		//}

		loaded_points++;
	}
	loaded++;
	fclose(f);
}
// ========================================================================================================================================
void PrepareVariables(){
	thetaarr  = new float*[5];
	etaarr    = new float*[5];
	phiarr    = new float*[5];
	x0arr     = new float*[5];
	x0arr_copy= new float*[5];
	x0Earr    = new float*[5];
	lamda0arr = new float*[5];

	for(int l = 0 ; l < loaded ; l++){
		thetaarr  [l] = new float[thetavec [l].size()];
		etaarr    [l] = new float[thetavec [l].size()];
		phiarr    [l] = new float[phivec   [l].size()];
		x0arr     [l] = new float[x0vec_max[l].size()];
		x0arr_copy[l] = new float[x0vec_max[l].size()];
		x0Earr    [l] = new float[x0vec_max[l].size()];
		lamda0arr [l] = new float[lamda0vec[l].size()];

		float x0min = 100.;

		for (int i=0; i<thetavec[l].size(); i++){
			thetaarr  [l][i] = 270. - thetavec[l][i];
			etaarr    [l][i] = -TMath::Log(TMath::Tan(thetaarr[l][i]/2.*TMath::Pi()/180.));
			phiarr    [l][i] = phivec[l][i];
			x0arr     [l][i] = (x0vec_max[l][i]+x0vec_min[l][i])/2.;
			x0Earr    [l][i] = (x0vec_max[l][i]-x0vec_min[l][i])/2.;
			x0arr_copy[l][i] = (x0vec_max[l][i]+x0vec_min[l][i])/2.;
			lamda0arr [l][i] = lamda0vec[l][i];
			if(etaarr[l][i]!=etaarr[l][i]){etaarr[l][i] = -10.; x0arr[l][i] = -10.;}
			if(x0arr [l][i] < x0min) x0min = x0arr[l][i];
		}

	}
}
// ========================================================================================================================================
void CreatePlots(){
	for(int l = 0 ; l < loaded ; l++){
		gr_th[l] = new TGraphErrors(thetavec[l].size(),thetaarr[l],x0arr_copy[l],0,x0Earr[l]);
		gr_et[l] = new TGraphErrors(thetavec[l].size(),etaarr  [l],x0arr     [l],0,x0Earr[l]);

		PrettyTGraphErrors( gr_th[l] , color[l] , 21 , -90, 270 , 999 , 999 , "#theta [deg]" , "radiation length (X/X_{0})" , "" );
		PrettyTGraphErrors( gr_et[l] , color[l] , 21 , -4 ,   4 ,   0 , .25 , "#eta"         , "radiation length (X/X_{0})" , "" );

		if(part[l]=="beampipe"){
			gr_th[l] -> SetLineStyle(2);
			gr_et[l] -> SetLineStyle(2);
		}
	}
}
// ========================================================================================================================================
void PrettyTGraphErrors( TGraphErrors * gP , int color , int marker , float xmin , float xmax , float ymin , float ymax , TString xtit , TString ytit , TString tit ){
	gP -> GetXaxis() -> SetRangeUser(xmin,xmax);
	gP -> SetMarkerSize(0.1);
	gP -> SetLineColor(color);
	gP -> SetFillColorAlpha(color,.3);
	gP -> SetLineWidth(3);
	gP -> SetMarkerColor(color);

	gP -> SetMinimum(0);

	if(ymax!=999) gP -> SetMaximum(ymax);
	if(ymin!=999) gP -> SetMinimum(ymin);

	gP -> SetTitle(tit);
	gP -> GetXaxis() -> SetTitle(xtit);
	gP -> GetYaxis() -> SetTitle(ytit);

	gP -> GetXaxis() -> CenterTitle();
	gP -> GetYaxis() -> CenterTitle();
}                
