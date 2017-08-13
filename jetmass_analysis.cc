//Esha Rao
//August 10, 2017
//Adding on even more to create jet shape plot based off of their jet mass to see for any correlation
//
//Esha Rao
//July 24, 2017
//Adding on to Aditya's code to create histograms of jet shape and jet mass
//
//
//Aditya Verma
//June 15, 2017
//This code creates a histigram of Z=pttrack/ptjet which is the fragmentation analysis
//
//
//
//
//
//

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class jetmass_analysis : public Analysis {
  public:

    /// Constructor
    jetmass_analysis()
      : Analysis("jetmass_analysis")
    {    }

    void init() {

      //! finalstate particles and the jets with the algorithm and radius definition
      FinalState fs(-5.0, 5.0, 0.150*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets");

   

      ChargedFinalState cfs(-5.0, 5.0, 0.100*GeV);
      addProjection(cfs, "CFS");

      //! initialize variables
      _ptCut = 10.0;
      _nJets = 0.0;
      _nevts = 0.0;
      _totJets=0;
      _twojets=0, _fourjets=0, _sixjets=0, _eightjets=0, _tenjets=0, _twelvejets=0, _fourteenjets=0;
     
      
      //! initialize histograms		      
      /* _h_Z = bookHisto1D("Z", 10, 0, 1);
	 _h_Aj = bookHisto1D("Aj", 20, 0, 1);
	 _h_delphi = bookHisto1D("Delta_Phi", 30, 1, 5);
	 _h_Iaa = bookHisto1D("Iaa",20, 0, 100);*/
      /* _h_jetshape1 = bookHisto1D("jetshape1", 20, 0, 1);
	 _h_jetshape2 = bookHisto1D("jetshape2", 20, 0, 1);
	 _h_jetshape3 = bookHisto1D("jetshape3", 20, 0, 1);
	 _h_jetshape4 = bookHisto1D("jetshape4", 20, 0, 1);
	 _h_jetshape5 = bookHisto1D("jetshape5", 20, 0, 1);
	 _h_jetshape6 = bookHisto1D("jetshape6", 20, 0, 1);
	 _h_jetshape7 = bookHisto1D("jetshape7", 20, 0, 1);
	 _h_jetshape8 = bookHisto1D("jetshape8", 20, 0, 1);*/
      _h_jetshape = bookHisto1D("jetshape", 8, 0, 0.4);
      _h_jetshape1 = bookHisto1D("jetshape1", 8, 0, 0.4);
      _h_jetshape2 = bookHisto1D("jetshape2", 8, 0, 0.4);
      _h_jetshape3 = bookHisto1D("jetshape3", 8, 0, 0.4);
      _h_jetshape4 = bookHisto1D("jetshape4", 8, 0, 0.4);
      _h_jetshape5 = bookHisto1D("jetshape5", 8, 0, 0.4);
      _h_jetshape6 = bookHisto1D("jetshape6", 8, 0, 0.4);
      _h_jetmass = bookHisto1D("jetmass", 8, 0, 80);    

      
    }

    //intialize global variables
    void analyze(const Event& event) {

      
      int nJetpE=0;
      double tracksum[8];

      double jetsum[8];

      jetsum[0]=0;
      jetsum[1]=0;
      jetsum[2]=0;
      jetsum[3]=0;
      jetsum[4]=0;
      jetsum[5]=0;
      jetsum[6]=0;
      jetsum[7]=0;

      double jetsum1[8];

      jetsum1[0]=0;
      jetsum1[1]=0;
      jetsum1[2]=0;
      jetsum1[3]=0;
      jetsum1[4]=0;
      jetsum1[5]=0;
      jetsum1[6]=0;
      jetsum1[7]=0;

      double jetsum2[8];

      jetsum2[0]=0;
      jetsum2[1]=0;
      jetsum2[2]=0;
      jetsum2[3]=0;
      jetsum2[4]=0;
      jetsum2[5]=0;
      jetsum2[6]=0;
      jetsum2[7]=0;

      double jetsum3[8];

      jetsum3[0]=0;
      jetsum3[1]=0;
      jetsum3[2]=0;
      jetsum3[3]=0;
      jetsum3[4]=0;
      jetsum3[5]=0;
      jetsum3[6]=0;
      jetsum3[7]=0;

      double jetsum4[8];

      jetsum4[0]=0;
      jetsum4[1]=0;
      jetsum4[2]=0;
      jetsum4[3]=0;
      jetsum4[4]=0;
      jetsum4[5]=0;
      jetsum4[6]=0;
      jetsum4[7]=0;      
	
      double jetsum5[8];

      jetsum5[0]=0;
      jetsum5[1]=0;
      jetsum5[2]=0;
      jetsum5[3]=0;
      jetsum5[4]=0;
      jetsum5[5]=0;
      jetsum5[6]=0;
      jetsum5[7]=0;

      double jetsum6[8];

      jetsum6[0]=0;
      jetsum6[1]=0;
      jetsum6[2]=0;
      jetsum6[3]=0;
      jetsum6[4]=0;
      jetsum6[5]=0;
      jetsum6[6]=0;
      jetsum6[7]=0;       
	
      //! event cross section weight 
      const double weight = event.weight(); 
      const ChargedFinalState tracks= applyProjection<ChargedFinalState>(event, "CFS");

      //! Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-1, 1) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);
      //   const Jet bjets = Ajets.jetsByPt(cuts);
      if(ajets.size()<2)
	vetoEvent;


      bool leadjetcut = false;
      if(ajets.at(0).pt() > 30)
	leadjetcut = true;

      bool dijetcut = false;
      if(ajets.at(0).pt() > 30 && ajets.at(1).pt() > 5){
	dijetcut = true;
	_nevts+=weight;
      
      }
      
      
      //! Fill histograms 
      //! Loop over jets in the event
      foreach ( PseudoJet jet, ajets ){
	
	_totJets++;

	tracksum[0]=0;
	tracksum[1]=0;
	tracksum[2]=0;
	tracksum[3]=0;
	tracksum[4]=0;
	tracksum[5]=0;
	tracksum[6]=0;
	tracksum[7]=0;
	
	int nTracks =0;
	_nJets+=weight;//Jet counter

	double delR;
	
	foreach ( Particle track, tracks.particles() ){
	    
	  delR = deltaR(track.eta(),track.phi(MINUSPI_PLUSPI),jet.eta(),jet.phi_std());	              

	  if(delR<0.05){
	    tracksum[0]+=track.pt()/jet.pt();}
	  else if(delR<0.1){
	    tracksum[1]+=track.pt()/jet.pt();}
	  else if(delR<0.15){
	    tracksum[2]+=track.pt()/jet.pt();}
          else if(delR<0.2){
	    tracksum[3]+=track.pt()/jet.pt();}
	  else if (delR<0.25){
	    tracksum[4]+=track.pt()/jet.pt();}
	  else if (delR<0.3){
	    tracksum[5]+=track.pt()/jet.pt();}
	  else if (delR<0.35){
	    tracksum[6]+=track.pt()/jet.pt();}
	  else if (delR<0.4){
	    tracksum[7]+=track.pt()/jet.pt();}	 

	}//! particle loop
	  
	//number of jets per event
	nJetpE++;

	//each different histogram for jet shape for a specific bin range from jet center
	/*      	_h_jetshape1->fill(tracksum[0], weight);
			_h_jetshape2->fill(tracksum[1], weight);
			_h_jetshape3->fill(tracksum[2], weight);
			_h_jetshape4->fill(tracksum[3], weight);
			_h_jetshape5->fill(tracksum[4], weight);
			_h_jetshape6->fill(tracksum[5], weight);
			_h_jetshape7->fill(tracksum[6], weight);
			_h_jetshape8->fill(tracksum[7], weight);*/

	//add up all the sum of track pt over jet pt in a specific region of jet mass  and put into the jet sum
	if(jet.m()<2){
	  jetsum[0]+=tracksum[0];
	  jetsum[1]+=tracksum[1];
	  jetsum[2]+=tracksum[2];		
	  jetsum[3]+=tracksum[3];
	  jetsum[4]+=tracksum[4];
	  jetsum[5]+=tracksum[5];
	  jetsum[6]+=tracksum[6];
	  jetsum[7]+=tracksum[7];
	  _twojets++;}
	else if(jet.m()>2 && jet.m()<4){
	  jetsum1[0]+=tracksum[0];
	  jetsum1[1]+=tracksum[1];
	  jetsum1[2]+=tracksum[2];		
	  jetsum1[3]+=tracksum[3];
	  jetsum1[4]+=tracksum[4];
	  jetsum1[5]+=tracksum[5];
	  jetsum1[6]+=tracksum[6];
	  jetsum1[7]+=tracksum[7];
	  _fourjets++;}
	else if(jet.m()>4 && jet.m()<6){
	  jetsum2[0]+=tracksum[0];
	  jetsum2[1]+=tracksum[1];
	  jetsum2[2]+=tracksum[2];		
	  jetsum2[3]+=tracksum[3];
	  jetsum2[4]+=tracksum[4];
	  jetsum2[5]+=tracksum[5];
	  jetsum2[6]+=tracksum[6];
	  jetsum2[7]+=tracksum[7];
	  _sixjets++;}
	else if(jet.m()>6 && jet.m()<8){
	  jetsum3[0]+=tracksum[0];
	  jetsum3[1]+=tracksum[1];
	  jetsum3[2]+=tracksum[2];		
	  jetsum3[3]+=tracksum[3];
	  jetsum3[4]+=tracksum[4];
	  jetsum3[5]+=tracksum[5];
	  jetsum3[6]+=tracksum[6];
	  jetsum3[7]+=tracksum[7];
	  _eightjets++;}
	else if(jet.m()>8 && jet.m()<10){
	  jetsum4[0]+=tracksum[0];
	  jetsum4[1]+=tracksum[1];
	  jetsum4[2]+=tracksum[2];		
	  jetsum4[3]+=tracksum[3];
	  jetsum4[4]+=tracksum[4];
	  jetsum4[5]+=tracksum[5];
	  jetsum4[6]+=tracksum[6];
	  jetsum4[7]+=tracksum[7];
	  _tenjets++;}
	else if(jet.m()>10 && jet.m()<12){
	  jetsum5[0]+=tracksum[0];
	  jetsum5[1]+=tracksum[1];
	  jetsum5[2]+=tracksum[2];		
	  jetsum5[3]+=tracksum[3];
	  jetsum5[4]+=tracksum[4];
	  jetsum5[5]+=tracksum[5];
	  jetsum5[6]+=tracksum[6];
	  jetsum5[7]+=tracksum[7];
	  _twelvejets++;}
	else if(jet.m()>12 && jet.m()<14){
	  jetsum6[0]+=tracksum[0];
	  jetsum6[1]+=tracksum[1];
	  jetsum6[2]+=tracksum[2];		
	  jetsum6[3]+=tracksum[3];
	  jetsum6[4]+=tracksum[4];
	  jetsum6[5]+=tracksum[5];
	  jetsum6[6]+=tracksum[6];
	  jetsum6[7]+=tracksum[7];
	  _fourteenjets++;}
	
      }//! jet loop

      //jet shape with bins of .05
      _h_jetshape->fill(0.025, jetsum[0]);
      _h_jetshape->fill(0.075, jetsum[1]);
      _h_jetshape->fill(0.125, jetsum[2]);
      _h_jetshape->fill(0.175, jetsum[3]);
      _h_jetshape->fill(0.225, jetsum[4]);
      _h_jetshape->fill(0.275, jetsum[5]);
      _h_jetshape->fill(0.325, jetsum[6]);
      _h_jetshape->fill(0.375, jetsum[7]);

      _h_jetshape1->fill(0.025, jetsum1[0]);
      _h_jetshape1->fill(0.075, jetsum1[1]);
      _h_jetshape1->fill(0.125, jetsum1[2]);
      _h_jetshape1->fill(0.175, jetsum1[3]);
      _h_jetshape1->fill(0.225, jetsum1[4]);
      _h_jetshape1->fill(0.275, jetsum1[5]);
      _h_jetshape1->fill(0.325, jetsum1[6]);
      _h_jetshape1->fill(0.375, jetsum1[7]);

      _h_jetshape2->fill(0.025, jetsum2[0]);
      _h_jetshape2->fill(0.075, jetsum2[1]);
      _h_jetshape2->fill(0.125, jetsum2[2]);
      _h_jetshape2->fill(0.175, jetsum2[3]);
      _h_jetshape2->fill(0.225, jetsum2[4]);
      _h_jetshape2->fill(0.275, jetsum2[5]);
      _h_jetshape2->fill(0.325, jetsum2[6]);
      _h_jetshape2->fill(0.375, jetsum2[7]);
	
      _h_jetshape3->fill(0.025, jetsum3[0]);
      _h_jetshape3->fill(0.075, jetsum3[1]);
      _h_jetshape3->fill(0.125, jetsum3[2]);
      _h_jetshape3->fill(0.175, jetsum3[3]);
      _h_jetshape3->fill(0.225, jetsum3[4]);
      _h_jetshape3->fill(0.275, jetsum3[5]);
      _h_jetshape3->fill(0.325, jetsum3[6]);
      _h_jetshape3->fill(0.375, jetsum3[7]);

      _h_jetshape4->fill(0.025, jetsum4[0]);
      _h_jetshape4->fill(0.075, jetsum4[1]);
      _h_jetshape4->fill(0.125, jetsum4[2]);
      _h_jetshape4->fill(0.175, jetsum4[3]);
      _h_jetshape4->fill(0.225, jetsum4[4]);
      _h_jetshape4->fill(0.275, jetsum4[5]);
      _h_jetshape4->fill(0.325, jetsum4[6]);
      _h_jetshape4->fill(0.375, jetsum4[7]);

      _h_jetshape5->fill(0.025, jetsum5[0]);
      _h_jetshape5->fill(0.075, jetsum5[1]);
      _h_jetshape5->fill(0.125, jetsum5[2]);
      _h_jetshape5->fill(0.175, jetsum5[3]);
      _h_jetshape5->fill(0.225, jetsum5[4]);
      _h_jetshape5->fill(0.275, jetsum5[5]);
      _h_jetshape5->fill(0.325, jetsum5[6]);
      _h_jetshape5->fill(0.375, jetsum5[7]);

      _h_jetshape6->fill(0.025, jetsum6[0]);
      _h_jetshape6->fill(0.075, jetsum6[1]);
      _h_jetshape6->fill(0.125, jetsum6[2]);
      _h_jetshape6->fill(0.175, jetsum6[3]);
      _h_jetshape6->fill(0.225, jetsum6[4]);
      _h_jetshape6->fill(0.275, jetsum6[5]);
      _h_jetshape6->fill(0.325, jetsum6[6]);
      _h_jetshape6->fill(0.375, jetsum6[7]);
      
    }//event loop

    /*  double shapeCalc(double R, double delR, double pT){

	if(dR<=R&&dR>=(R-0.1))
	return pT;
	else
	return 0.0;

      
      
	}
    */
    /*   void fillThis( double a[]){

	 _h_jetshape->fill(0.05, a[0]/_toJets);
	 _h_jetshape->fill(0.15, a[1]/_totJets);
	 _h_jetshape->fill(0.25, a[2]/_totJets);
	 _h_jetshape->fill(0.35, a[3]/_totJets);

	 }
    */

    /// Normalise histograms etc., after the run
    void finalize() {/*


_h_jetshape->fill(0.05, tracksum[0]/jet.pt());
_h_jetshape->fill(0.15, tracksum[1]/jet.pt());
_h_jetshape->fill(0.25, tracksum[2]/jet.pt());
_h_jetshape->fill(0.35, tracksum[3]/jet.pt());
		     */
      
      /*    scale(_h_jetpT, crossSection()/picobarn/sumOfWeights());
	    scale(_h_Z, 1./_nJets);
	    scale(_h_Aj, 1./_nevts);
	    scale(_h_delphi, 1./_nevts);
	    scale(_h_Iaa,crossSection()/picobarn/sumOfWeights());*/
      scale(_h_jetmass, 1./_totJets);
      scale(_h_jetshape, 1./_twojets);
      scale(_h_jetshape1, 1./_fourjets);
      scale(_h_jetshape2, 1./_sixjets);
      scale(_h_jetshape3, 1./_eightjets);
      scale(_h_jetshape4, 1./_tenjets);
      scale(_h_jetshape5, 1./_twelvejets);
      scale(_h_jetshape6, 1./_fourteenjets);
    }

    

    //! Declare variables 
  private:
    /* Histo1DPtr _h_jetshape1;
       Histo1DPtr _h_jetshape2;
       Histo1DPtr _h_jetshape3;
       Histo1DPtr _h_jetshape4;
       Histo1DPtr _h_jetshape5;
       Histo1DPtr _h_jetshape6;
       Histo1DPtr _h_jetshape7;
       Histo1DPtr _h_jetshape8;*/
    Histo1DPtr _h_jetshape1;
    Histo1DPtr _h_jetshape2;
    Histo1DPtr _h_jetshape3;
    Histo1DPtr _h_jetshape4;
    Histo1DPtr _h_jetshape5;
    Histo1DPtr _h_jetshape6;
    Histo1DPtr _h_jetshape;
    Histo1DPtr _h_jetmass;

  protected:
    double _ptCut;
    double _nJets;//Jet Counter, Track Counter
    double _nevts;
    double _totJets;
    int _twojets, _fourjets, _sixjets, _eightjets, _tenjets, _twelvejets, _fourteenjets; 
    // double _trackSum[_totJets];
    
  };


  // The hook for the plugin syst
  DECLARE_RIVET_PLUGIN(jetmass_analysis);


}
