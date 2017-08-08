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
      _ptCut = 30.0;
      _nJets = 0.0;
      _nevts = 0.0;
      _totJets=0;
     
      
      //! initialize histograms		      
      /* _h_Z = bookHisto1D("Z", 10, 0, 1);
      _h_Aj = bookHisto1D("Aj", 20, 0, 1);
      _h_delphi = bookHisto1D("Delta_Phi", 30, 1, 5);
      _h_Iaa = bookHisto1D("Iaa",20, 0, 100);*/
      _h_jetshape1 = bookHisto1D("jetshape1", 20, 0, 1);
      _h_jetshape2 = bookHisto1D("jetshape2", 20, 0, 1);
      _h_jetshape3 = bookHisto1D("jetshape3", 20, 0, 1);
      _h_jetshape4 = bookHisto1D("jetshape4", 20, 0, 1);
      _h_jetshape5 = bookHisto1D("jetshape5", 20, 0, 1);
      _h_jetshape6 = bookHisto1D("jetshape6", 20, 0, 1);
      _h_jetshape7 = bookHisto1D("jetshape7", 20, 0, 1);
      _h_jetshape8 = bookHisto1D("jetshape8", 20, 0, 1);
      _h_jetshape = bookHisto1D("jetshape", 8, 0, 0.4);
      _h_jetmass = bookHisto1D("jetmass", 8, 0, 80)    

      
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

	//	_h_jetmass->fill(bjets.mass(), weight);
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

	  /*  tracksum[0]+=shapeCalc(0.1, delR, track.pt())/jet.pt();
	  tracksum[1]+=shapeCalc(0.2, delR, track.pt())/jet.pt();
	  tracksum[2]+=shapeCalc(0.3, delR, track.pt())/jet.pt();
	  tracksum[3]+=shapeCalc(0.4, delR, track.pt())/jet.pt();
	  */
	}//! particle loop

	//Jet Pt
	//_h_jetpT->fill(jet.pt(), weight);

	//Jet Eta
	//	_h_jeteta->fill(jet.eta(), weight);

	//Jet Phi
	//_h_jetphi->fill(jet.phi_std(), weight);
	  
	
	nJetpE++;

      	_h_jetshape1->fill(tracksum[0], weight);
	_h_jetshape2->fill(tracksum[1], weight);
	_h_jetshape3->fill(tracksum[2], weight);
	_h_jetshape4->fill(tracksum[3], weight);
	_h_jetshape5->fill(tracksum[4], weight);
	_h_jetshape6->fill(tracksum[5], weight);
	_h_jetshape7->fill(tracksum[6], weight);
	_h_jetshape8->fill(tracksum[7], weight);

       	jetsum[0]+=tracksum[0];
	jetsum[1]+=tracksum[1];
	jetsum[2]+=tracksum[2];		
	jetsum[3]+=tracksum[3];
	jetsum[4]+=tracksum[4];
	jetsum[5]+=tracksum[5];
	jetsum[6]+=tracksum[6];
	jetsum[7]+=tracksum[7];

      //seperate mass plots
      if(jet.m()>30 && jet.m()<35){
      }
	
      }//! jet loop
    
    	_h_jetshape->fill(0.025, jetsum[0]);
	_h_jetshape->fill(0.075, jetsum[1]);
        _h_jetshape->fill(0.125, jetsum[2]);
	_h_jetshape->fill(0.175, jetsum[3]);
	_h_jetshape->fill(0.225, jetsum[4]);
	_h_jetshape->fill(0.275, jetsum[5]);
	_h_jetshape->fill(0.325, jetsum[6]);
	_h_jetshape->fill(0.375, jetsum[7]);


      
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
      scale(_h_jetshape, 1./_totJets);
    }

    

    //! Declare variables 
  private:
    Histo1DPtr _h_jetshape1;
    Histo1DPtr _h_jetshape2;
    Histo1DPtr _h_jetshape3;
    Histo1DPtr _h_jetshape4;
    Histo1DPtr _h_jetshape5;
    Histo1DPtr _h_jetshape6;
    Histo1DPtr _h_jetshape7;
    Histo1DPtr _h_jetshape8;
    Histo1DPtr _h_jetshape;
    Histo1DPtr _h_jetmass;

  protected:
    double _ptCut;
    double _nJets;//Jet Counter, Track Counter
    double _nevts;
    double _totJets;
    // double _trackSum[_totJets];
    
  };


  // The hook for the plugin syst
  DECLARE_RIVET_PLUGIN(jetmass_analysis);


}
