//Esha Rao
//August 2, 2017
//Writing a code from scratch for basic properties of jets such as pT, phi, eta (will add cuts later based on plots)

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class jetproperties_analysis : public Analysis {
  public:

    /// Constructor
    jetproperties_analysis()
      : Analysis("jetproperties_analysis")
    {    }

//----------------------------------------------------------------------------------------------------------------
    
    // Book histograms and initialise projections before the run
    void init() {

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

      //! initialize histograms
      _h_jetpT = bookHisto1D("jetpT", 20, 0, 100);
      _h_jeteta = bookHisto1D("jeteta", 100, -2.5, 2.5);
      _h_jetphi = bookHisto1D("jetphi", 100, -3.5, 3.5);
      _h_jetmass = bookHisto1D("jetmass",50, 0, 100);
      _h_jetphivseta = bookHisto2D("jetphivseta",100, -2.5, 2.5, 100, -3.5, 3.5);
      
    }

//-----------------------------------------------------------------------------------------------------------------   

    // Perform the per-event analysis
    void analyze(const Event& event) {

      //Number of jets per event
      int nJetpE=0;
      
      //Array for sum of pt of track divided pt of jet for a specific range from jet axis
      // double tracksum[8];

      //! event cross section weight 
      const double weight = event.weight(); 
      const ChargedFinalState tracks= applyProjection<ChargedFinalState>(event, "CFS");

      // Jet collection of importance 
      const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      Cut cuts = Cuts::etaIn(-1, 1) & (Cuts::pT > _ptCut*GeV);
      const Jets ajets = Ajets.jetsByPt(cuts);

    //Jet loop
    foreach(PseudoJet jet, ajets){

      //Total Number of Jets
      _totJets++;

      //Every time code loops over a jet again the array will be cleared
      /*tracksum[0]=0;
      tracksum[1]=0;
      tracksum[2]=0;
      tracksum[3]=0;
      tracksum[4]=0;
      tracksum[5]=0
      tracksum[6]=0;
      tracksum[7]=0;*/

      double delR;
      int nTracks=0;

      //Jet counter
      _nJets+=weight;

      //Track Loop
      foreach(Particle track, tracks.particles()){

	//Delta R equals
	delR = deltaR(track.eta(),track.phi(MINUSPI_PLUSPI),jet.eta(),jet.phi_std());
	
      }//End of Track Loop

      //Jet Pt
      _h_jetpT->fill(jet.pt(), weight);

      //Jet Eta
      _h_jeteta->fill(jet.eta(), weight);

      //Jet Phi
      _h_jetphi->fill(jet.phi_std(), weight);

      //Jet Mass
      _h_jetmass->fill(jet.m(), weight);

      //Jet Phi vs Eta
      _h_jetphivseta->fill(jet.eta(), jet.phi(), weight);
      
    }//End of Jet Loop

   }//End of Event Loop
//----------------------------------------------------------------------------------------------------------------

    //Normalise histograms etc., after the run
    void finalize(){
      scale(_h_jetmass, 1./_nJets);
    }

//----------------------------------------------------------------------------------------------------------------

    //Declare variables
  private:
    Histo1DPtr _h_jetpT;
    Histo1DPtr _h_jeteta;
    Histo1DPtr _h_jetphi;
    Histo1DPtr _h_jetmass;
    Histo2DPtr _h_jetphivseta;
    
  protected:
    double _ptCut;
    double _nJets;//Jet Counter
    double _nevts;
    double _totJets;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(jetproperties_analysis);

}

