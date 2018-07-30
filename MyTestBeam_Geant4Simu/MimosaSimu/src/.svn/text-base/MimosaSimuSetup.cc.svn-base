// @(#)maf/dtools:$Name:  $:$Id:MimosaSimuSetup.cxx  v.1 2005/10/02 18:03:46 sha Exp $

#include "MimosaSimuSetup.hh"
#include "TApplication.h"
#include "TFile.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TVector2.h"

#include <assert.h>

//ClassImp(MimosaSimuSetup) // MimosaSimuSetup of Test Beam Setup and Data Structure

//==================================================================
MimosaSimuSetup::MimosaSimuSetup(){ 
  MimosaSimuSetupDebug = 0;
  fFieldMaxLength = 100;
  fFieldName = new Char_t[fFieldMaxLength];
  
  InitSourceList();
}
//==================================================================  
MimosaSimuSetup::MimosaSimuSetup(const char* ConfigFile)
{ 

  // Modified JB 2012/12/20 new member data "fieldName"

  MimosaSimuSetupDebug = 0;
  fConfigPathAndFileName = TString(ConfigFile);
  fFieldMaxLength = 100;
  fFieldName = new Char_t[fFieldMaxLength];

  InitSourceList();

}
//==================================================================
void MimosaSimuSetup::SetConfigFileName(TString aCFN) 
{ 
  fConfigPathAndFileName = aCFN;     
  if(MimosaSimuSetupDebug) cout << "MimosaSimuSetup::SetConfigFileName "<< aCFN << endl ;
  
  InitSourceList();
  
} 
//==================================================================
void MimosaSimuSetup::ReadRunParameters() 
{  
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Run Parameter 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified JB 2013/08/19: default parameters setup
  // Modified VR 2014/06/30: RunNumber is not read anymore, taken from DSession 
  // Modified VR 2014/06/30: Add the DataSubDirPrefix and test of mandatory config
  
  cout << endl << " - Reading Run parameters" << endl;
  
  TString RunParInitString("");
  // Set default values for optional or mandatory parameters
  RunParameter.Affiliation      = RunParInitString;
  RunParameter.Signature        = RunParInitString;
  RunParameter.BeamTime         = RunParInitString;
  RunParameter.Confidence       = RunParInitString;
  RunParameter.DataPath         = TString("-1");// [MANDATORY]
  RunParameter.DataSubDirPrefix = RunParInitString;
  RunParameter.Extension        = RunParInitString;
  
  RunParameter.Number           = 0;
  RunParameter.EventsInFile     = 0;// [MANDATORY for IMGBoardReader]
  RunParameter.StartIndex       = 0;// [MANDATORY for IMGBoardReader]
  RunParameter.EndIndex         = 0;// [MANDATORY for IMGBoardReader]
  RunParameter.NoiseRun         = 0;
  
  do {  
    nextField();

    if ( ! strcmp( fFieldName, "Affiliation" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.Affiliation = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "Signature" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.Signature = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "BeamTime" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.BeamTime = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "Confidence" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.Confidence = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "DataPath" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.DataPath = TString(mmm);
      //sprintf( RunParameter.DataPath, "%s", "");
      RunParameter.DataPath = TString("");
    }
    else if( ! strcmp( fFieldName, "DataSubDirPrefix" ) ) { // VR 2014/06/30 
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.DataSubDirPrefix = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "Extension" ) ) {
      Char_t mmm[RunParameter.tpsz];
      read_strings( mmm, RunParameter.tpsz);
      RunParameter.Extension = TString(mmm);
    }
    //VR 2014/06/30 RunNumber is set from DSession
    else if( ! strcmp( fFieldName, "RunNumber" ) ) {
      cout << "WARNING : parameter 'RunNumber' in config file is ignored, because it is redundant with InitSession command argument !" << endl;
      getRidOfLine();
    }
    else if( ! strcmp( fFieldName, "EventsInFile" ) ) {
      read_item(RunParameter.EventsInFile);
    }
    else if( ! strcmp( fFieldName, "StartIndex" ) ) {
      read_item(RunParameter.StartIndex);
    }
    else if( ! strcmp( fFieldName, "EndIndex" ) ) {
      read_item(RunParameter.EndIndex);
    }
    else if( ! strcmp( fFieldName, "NoiseRun" ) ) {
      read_item(RunParameter.NoiseRun);       //YV 27/11/09
    }
    else
    {
      if ( strcmp( fFieldName, "Planes") && strcmp( fFieldName, "Ladders") )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }
  } while ( strcmp( fFieldName, "Planes") && strcmp( fFieldName, "Ladders") );
  
#if 0
  // Test of MANDATORY field settings
  // VR 20014/06/30 
  if ( !strcmp(RunParameter.DataPath,"-1"))
  {
    cout << " ERROR: The field 'DataPath' is MANDATORY in config file (section run parameters)" << endl;
    gApplication->Terminate();
  }
#endif
}
//==================================================================
void MimosaSimuSetup::ReadTrackerParameters() 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Tracker 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified: LC 2013/??/?? new parameters for vertex finder
  // Modified: JB 2013/06/21 new parameters to control track_finder
  // Modified: JB 2013/07/17 new parameters to select track_finder method
  // Modified: VR 2014/06/29 add the SearchMoreHitDistance parameter and test of mandatory config
  // Modified: VR 2014/10/30 add the PreTrackHitsNbMinimum and HitsNbMinForFullTrack parameters and test of mandatory config
  // Modified: VR 2014/12/18 modify all parameters for track_finder 2 (IVI) 
  // Modified: JB 2014/12/23 add Subtrack parameters
  // Modified: LC 2015/01/?? add HitMonteCarlo
  // Modified: LC 2015/01/31 add DprecAlignMethod
   
  cout << endl << " - Reading Parameters of the Tracker" << endl;

  // ################################################
  //            Set default values 
  // ################################################
  // *****************************
  //       Planes config
  // *****************************
  TrackerParameter.Planes             = 0;
  TrackerParameter.Ladders            = 0;
  TrackerParameter.TimeLimit          = -1;
  TrackerParameter.Resolution         = -1.;
  TrackerParameter.ExpSetup           = TString("");  //Default experimental setup is beam-test
  //Beam-Test experimental setup parameters
  TrackerParameter.BeamNparticles     = -1.0;
  TrackerParameter.BeamRandNparticles = false;
  TrackerParameter.BeamType           = TString("");
  TrackerParameter.BeamMomentum       = -1.0;
  TrackerParameter.BeamDirection      = G4ThreeVector(0.0,0.0,1.0);
  TrackerParameter.BeamAngularSpread  = G4ThreeVector(0.0,0.0,1.0);
  TrackerParameter.BeamMomentumSpread = -1.0;
  TrackerParameter.BeamOrigin         = G4ThreeVector(0.0,0.0,-1.0e+6);
  TrackerParameter.BeamOriginSpread   = G4ThreeVector(-1.0,-1.0,0.0);
  
  //Source experimental setup parameters
  TrackerParameter.SourcePosition     = G4ThreeVector(-999.0,-999.0,-999.0);
  TrackerParameter.SourceTilt         = G4ThreeVector(-999.0,-999.0,-999.0);
  TrackerParameter.SourceRadius       = -1.0;
  TrackerParameter.SourceActivity     = -1.0;
  TrackerParameter.SourceSensorROTime = -1.0;
  
  TrackerParameter.FillNonSensitiveBranch = false;
  
  //Magnetic field, assumed constant all over the world volume
  TrackerParameter.BFieldMagnitude    = 0.0;
  TrackerParameter.BFieldDirection    = G4ThreeVector(0.0,0.0,1.0);
  
  TrackerParameter.MediumMaterial = TString("");
  TrackerParameter.HitsInPlaneMaximum = -1 ;//[MANDATORY]
  
  TrackerParameter.TracksMaximum = 0;
  TrackerParameter.TracksFinder = 0; 
  TrackerParameter.KeepUnTrackedHitsBetw2evts = 0 ;
  TrackerParameter.PlanesForTrackMinimum = -1;//[MANDATORY if TracksMaximum!=0 && TracksFinder==0 ]
  TrackerParameter.HitsInPlaneTrackMaximum = -1;//[MANDATORY if TracksMaximum!=0 && TracksFinder==0 ]
  TrackerParameter.HitsInPlaneMaximum = -1;
  TrackerParameter.SearchHitDistance = -1.;//[MANDATORY if TracksMaximum<2] 
  TrackerParameter.SearchMoreHitDistance = -1.;//[MANDATORY if TracksMaximum<2]
  TrackerParameter.UseSlopeInTrackFinder = 1;
  TrackerParameter.TrackingPlaneOrderType = 0;
  TrackerParameter.SubtrackNplanes = 0; // JB 2014/12/23
  TrackerParameter.SubtrackPlanes = NULL;
  TrackerParameter.HitMonteCarlo = 0; // LC 2015/01
  TrackerParameter.DPrecAlignMethod = 0; // LC 2015/01/31
  
  TrackerParameter.Resolution     = 4.0; //um
  //TrackerParameter.BeamType       = TString("");
  //TrackerParameter.BeamMomentum   = 0.0; //GeV
  //TrackerParameter.MediumMaterial = TString("");
  
  // *****************************
  //  Tracking with track_finder 2
  // *****************************
  // all MANDATORY
  TrackerParameter.TrackingPass = -1;
  TrackerParameter.PreTrackHitsNbMinimum  = NULL;
  TrackerParameter.PreTrackHitsTypeUsed   = NULL;
  TrackerParameter.PreTrackHitsMaxDist    = NULL;
  TrackerParameter.ExtTrackHitsNbMinimum  = NULL;
  TrackerParameter.ExtTrackHitsTypeUsed   = NULL;
  TrackerParameter.ExtTrackHitsMaxDist    = NULL; 
  TrackerParameter.FullTrackHitsNbMinimum = NULL;
  
  TrackerParameter.TrackingOrderLines     = -1;
  TrackerParameter.TrackingOrderPreTrack = NULL;
  TrackerParameter.TrackingOrderExtTrack = NULL;
  
  // *****************************
  //          Vertexing
  // *****************************
  TrackerParameter.VertexMaximum = 0; // no vertexing by default
  TrackerParameter.VertexConstraint = 0; // use a vertex constraint to start track
  // *****************************
  //         Alignment
  // *****************************
  TrackerParameter.EventsForAlignmentFit = 0;
  
  do {
    // ----------------------------------
    //     Planes and their parameters
    // ----------------------------------
    if ( ! strcmp( fFieldName, "Planes") ) {
      read_item(TrackerParameter.Planes);  
    }
    else if ( ! strcmp( fFieldName, "Ladders") ) {
      read_item(TrackerParameter.Ladders);  
    }
    else if( ! strcmp( fFieldName, "TimeLimit" ) ) {
      read_item(TrackerParameter.TimeLimit); //RDM260609  
    }
    else if( ! strcmp( fFieldName, "Resolution" ) ) {
      read_item(TrackerParameter.Resolution); // RDM250809  
    }
    else if( ! strcmp( fFieldName, "ExpSetup" ) ) {  // AP 2016/07/08, Experimental setup: e.g. Beam-Test, sources, other ...
      read_TStrings( TrackerParameter.ExpSetup, 350);
    }
    //Beam-Test experimental setup parameters
    else if( ! strcmp( fFieldName, "BeamNparticles" ) ) {
      read_item(TrackerParameter.BeamNparticles);
    }
    else if( ! strcmp( fFieldName, "BeamRandNparticles" ) ) {
      int temp;
      read_item(temp);
      if(temp == 1) TrackerParameter.BeamRandNparticles = true;
      else          TrackerParameter.BeamRandNparticles = false;
    }
    else if( ! strcmp( fFieldName, "BeamMomentum" ) ) {
      read_item(TrackerParameter.BeamMomentum); // 2015/03/10, AP Beam momentum in GeV/c
    }
    else if( ! strcmp( fFieldName, "BeamType" ) ) {  // 2015/03/10, AP Beam type: e.g. e+/-, pi+/-, ...
      read_TStrings( TrackerParameter.BeamType, 350);
    }
    else if( ! strcmp( fFieldName, "BeamDirectionX" ) ) {
      read_r3(TrackerParameter.BeamDirection);
      if(MimosaSimuSetupDebug) cout << "  BeamDirection X = " << TrackerParameter.BeamDirection(0) << ",  Y = " << TrackerParameter.BeamDirection(1) << ", Z = " << TrackerParameter.BeamDirection(1) << endl;
    }
    else if( ! strcmp( fFieldName, "BeamAngularSpreadX" ) ) {
      Float_t temporal;
      read_item(temporal);
      TrackerParameter.BeamAngularSpread.setX(temporal);
    }
    else if( ! strcmp( fFieldName, "BeamAngularSpreadY" ) ) {
      Float_t temporal;
      read_item(temporal);
      TrackerParameter.BeamAngularSpread.setY(temporal);
    }
    else if( ! strcmp( fFieldName, "BeamMomentumSpread" ) ) {
      read_item(TrackerParameter.BeamMomentumSpread);
    }
    else if( ! strcmp( fFieldName, "BeamOriginX" ) ) {
      read_r3(TrackerParameter.BeamOrigin);
      TrackerParameter.BeamOrigin *= 1000; // move from millimeters to microns
      if(MimosaSimuSetupDebug) cout << "  BeamOrigin X = " << TrackerParameter.BeamOrigin(0) << ",  Y = " << TrackerParameter.BeamOrigin(1) << ", Z = " << TrackerParameter.BeamOrigin(1) << endl;
    }
    else if( ! strcmp( fFieldName, "BeamOriginSpreadX" ) ) {
      Float_t temporal;
      read_item(temporal);
      TrackerParameter.BeamOriginSpread.setX(temporal*1000);  // move from millimeters to microns
    }
    else if( ! strcmp( fFieldName, "BeamOriginSpreadY" ) ) {
      Float_t temporal;
      read_item(temporal);
      TrackerParameter.BeamOriginSpread.setY(temporal*1000);  // move from millimeters to microns
    }
    //Source experimental setup parameters
    else if( ! strcmp( fFieldName, "SourcePositionX" ) ) {
      read_r3(TrackerParameter.SourcePosition);
      TrackerParameter.SourcePosition *= 1000.0;  // move from millimeters to microns
      if(MimosaSimuSetupDebug) cout << "  SourcePosition X = " << TrackerParameter.SourcePosition(0) << ",  Y = " << TrackerParameter.SourcePosition(1) << ", Z = " << TrackerParameter.SourcePosition(1) << endl;
    }
    else if( ! strcmp( fFieldName, "SourceTiltZ" ) ) {
      read_r3(TrackerParameter.SourceTilt);
      TrackerParameter.SourceTilt *= M_PI/180.; // conversion degree to radian
      if(MimosaSimuSetupDebug) cout << "  SourceTilt X = " << TrackerParameter.SourcePosition(2) << ",  Y = " << TrackerParameter.SourceTilt(1) << ", Z = " << TrackerParameter.SourceTilt(0) << endl;
    }
    else if( ! strcmp( fFieldName, "SourceRadius" ) ) {
      read_item(TrackerParameter.SourceRadius);
      TrackerParameter.SourceRadius *= 1000.0; // move from millimeters to microns
    }
    else if( ! strcmp( fFieldName, "SourceActivity" ) ) {
      read_item(TrackerParameter.SourceActivity);
    }
    else if( ! strcmp( fFieldName, "SourceSensorROTime" ) ) {
      read_item(TrackerParameter.SourceSensorROTime);
    }
    else if( ! strcmp( fFieldName, "FillNonSensitiveBranch" ) ) {
      int temp_int;
      read_item(temp_int);
      if(temp_int == 1) TrackerParameter.FillNonSensitiveBranch = true;
    }
    //Magnetic field, assumed constant all over the world volume
    else if( ! strcmp( fFieldName, "BFieldMagnitude" ) ) {
      read_item(TrackerParameter.BFieldMagnitude); // 2016/07/08, AP B-field magnitude
    }
    else if( ! strcmp( fFieldName, "BFieldDirectionX" ) ) {
      read_r3(TrackerParameter.BFieldDirection);
      if(MimosaSimuSetupDebug) cout << "  BField direction X = " << TrackerParameter.BFieldDirection(0) << ",  Y = " << TrackerParameter.BFieldDirection(1) << ", Z = " << TrackerParameter.BFieldDirection(1) << endl;
    }
    
    else if( ! strcmp( fFieldName, "MediumMaterial" ) ) {  // 2015/06/01, AP Material of the medium containing the sensors: e.g. Dry air, or Vacuum
      read_TStrings( TrackerParameter.MediumMaterial, 350);
    }
    
    // ----------------------------------
    //     Clustering in planes
    // ----------------------------------
    else if( ! strcmp( fFieldName, "HitsInPlaneMaximum" ) ) {
      read_item(TrackerParameter.HitsInPlaneMaximum);  
    }
    // ----------------------------------
    //     Tracking parameters
    // ----------------------------------
    else if( ! strcmp( fFieldName, "TracksMaximum" ) ) {
      read_item(TrackerParameter.TracksMaximum);  
    }
    else if( ! strcmp( fFieldName, "TracksFinder" ) ) {
      read_item(TrackerParameter.TracksFinder); 
      cout << "    -> TracksFinder  = " << TrackerParameter.TracksFinder << endl ;
    }
    else if( ! strcmp( fFieldName, "PlanesForTrackMinimum" ) ) {
      read_item(TrackerParameter.PlanesForTrackMinimum);  
    }    
    else if( ! strcmp( fFieldName, "HitsInPlaneTrackMaximum" ) ) {
      read_item(TrackerParameter.HitsInPlaneTrackMaximum);
    }
    else if( ! strcmp( fFieldName, "SearchHitDistance" ) ) {
      read_item(TrackerParameter.SearchHitDistance);    
    }   
    else if( ! strcmp( fFieldName, "SearchMoreHitDistance" ) ) {
      read_item(TrackerParameter.SearchMoreHitDistance);    
    }
    else if( ! strcmp( fFieldName, "UseSlopeInTrackFinder" ) ) { // JB 2013/06/21
      read_item(TrackerParameter.UseSlopeInTrackFinder);  
    }
    else if( ! strcmp( fFieldName, "TrackingPlaneOrderType" ) ) { // VR 2014/07/14
      read_item(TrackerParameter.TrackingPlaneOrderType);  
    }
    else if( ! strcmp( fFieldName, "KeepUnTrackedHitsBetw2evts" ) ) { // VR 2014/08/26
      read_item(TrackerParameter.KeepUnTrackedHitsBetw2evts);  
    }
    else if( ! strcmp( fFieldName, "HitMonteCarlo" ) ) {
      read_item(TrackerParameter.HitMonteCarlo);    
    }
    else if( ! strcmp( fFieldName, "DPrecAlignMethod" ) ) {
      read_item(TrackerParameter.DPrecAlignMethod);    
    }
 
    // -------------------------------------------
    //     Tracking parameters for track_finder 2
    // -------------------------------------------
    // -*-*-*- Tracking Pass -*-*-*-
    else if( ! strcmp( fFieldName, "TrackingPass" ) )
    {
      read_item(TrackerParameter.TrackingPass);
      cout << "    -> Tracking pass = " << TrackerParameter.TrackingPass << endl ;

      if(TrackerParameter.TrackingPass>0)
      {
        Int_t nbOfReadingLoop = 0 ;
        while(kTRUE)
        {
          // Exit the loop
          if(TrackerParameter.PreTrackHitsNbMinimum  && \
            TrackerParameter.PreTrackHitsTypeUsed   && \
            TrackerParameter.PreTrackHitsMaxDist    && \
            TrackerParameter.ExtTrackHitsNbMinimum  && \
            TrackerParameter.ExtTrackHitsTypeUsed   && \
            TrackerParameter.ExtTrackHitsMaxDist    && \
            TrackerParameter.FullTrackHitsNbMinimum)
          {
            cout << "      -> Reading finished ! " << endl ;
            break;
          }
          
          if(nbOfReadingLoop > 20)
          {
            cout << "      -> ERROR: There is (a) missing parameter(s) : " << endl;
            if(! TrackerParameter.PreTrackHitsNbMinimum)  cout << "          * PreTrackHitsNbMinimum  list is not defined !" << endl ;
            if(! TrackerParameter.PreTrackHitsTypeUsed)   cout << "          * PreTrackHitsTypeUsed   list is not defined !" << endl ;
            if(! TrackerParameter.PreTrackHitsMaxDist)    cout << "          * PreTrackHitsMaxDist    list is not defined !" << endl ;
            if(! TrackerParameter.ExtTrackHitsNbMinimum)  cout << "          * ExtTrackHitsNbMinimum  list is not defined !" << endl ;
            if(! TrackerParameter.ExtTrackHitsTypeUsed)   cout << "          * ExtTrackHitsTypeUsed   list is not defined !" << endl ;
            if(! TrackerParameter.ExtTrackHitsMaxDist)    cout << "          * ExtTrackHitsMaxDist    list is not defined !" << endl ;
            if(! TrackerParameter.FullTrackHitsNbMinimum) cout << "          * FullTrackHitsNbMinimum list is not defined !" << endl ;
          
            gApplication->Terminate();
          }

          nextField();
          nbOfReadingLoop++;
	  
          if ( ! strcmp( fFieldName, "PreTrackHitsNbMinimum") ) 
          {
            TrackerParameter.PreTrackHitsNbMinimum  = new Int_t[TrackerParameter.TrackingPass];
            cout << "      -> PreTrackHitsNbMinimum {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.PreTrackHitsNbMinimum[iPass-1]);
              cout << " " << TrackerParameter.PreTrackHitsNbMinimum[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "PreTrackHitsTypeUsed") ) 
          {
            TrackerParameter.PreTrackHitsTypeUsed   = new Int_t[TrackerParameter.TrackingPass];
            cout << "      -> PreTrackHitsTypeUsed {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.PreTrackHitsTypeUsed[iPass-1]);
              cout << " " << TrackerParameter.PreTrackHitsTypeUsed[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "PreTrackHitsMaxDist") ) 
          {
            TrackerParameter.PreTrackHitsMaxDist    = new Float_t[TrackerParameter.TrackingPass];
            cout << "      -> PreTrackHitsMaxDist {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.PreTrackHitsMaxDist[iPass-1]);
              cout << " " << TrackerParameter.PreTrackHitsMaxDist[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "ExtTrackHitsNbMinimum") ) 
          {
            TrackerParameter.ExtTrackHitsNbMinimum  = new Int_t[TrackerParameter.TrackingPass];
            cout << "      -> ExtTrackHitsNbMinimum {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.ExtTrackHitsNbMinimum[iPass-1]);
              cout << " " << TrackerParameter.ExtTrackHitsNbMinimum[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "ExtTrackHitsTypeUsed") ) 
          {
            TrackerParameter.ExtTrackHitsTypeUsed   = new Int_t[TrackerParameter.TrackingPass];
            cout << "      -> ExtTrackHitsTypeUsed {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.ExtTrackHitsTypeUsed[iPass-1]);
              cout << " " << TrackerParameter.ExtTrackHitsTypeUsed[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "ExtTrackHitsMaxDist") ) 
          {
            TrackerParameter.ExtTrackHitsMaxDist    = new Float_t[TrackerParameter.TrackingPass];
            cout << "      -> ExtTrackHitsMaxDist {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.ExtTrackHitsMaxDist[iPass-1]);
              cout << " " << TrackerParameter.ExtTrackHitsMaxDist[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
          else if ( ! strcmp( fFieldName, "FullTrackHitsNbMinimum") ) 
          {
            TrackerParameter.FullTrackHitsNbMinimum = new Int_t[TrackerParameter.TrackingPass];
            cout << "      -> FullTrackHitsNbMinimum {" ; 
            nextItem('{');
            for(Int_t iPass=1 ; iPass <=TrackerParameter.TrackingPass ; iPass++)
            {
              read_item(TrackerParameter.FullTrackHitsNbMinimum[iPass-1]);
              cout << " " << TrackerParameter.FullTrackHitsNbMinimum[iPass-1] << " ";
            }
            cout << "}" << endl;
            nextItem('}');
          }
        }
      }
      else
      {
        cout << "      -> ERROR: The field 'TrackingPass' must be >0 ! " << endl;
        gApplication->Terminate();
      }    
    } // end if tracking pass
    
    
    // -*-*-*- Tracking Order -*-*-*-
    else if( ! strcmp( fFieldName, "TrackingOrderLines" ) )
    {
      read_item(TrackerParameter.TrackingOrderLines);
      cout << "    -> TrackingOrderLines = " << TrackerParameter.TrackingOrderLines << endl ;
      if(TrackerParameter.TrackingOrderLines>0)
      {
        TrackerParameter.TrackingOrderPreTrack = new Int_t*[TrackerParameter.TrackingOrderLines];
        TrackerParameter.TrackingOrderExtTrack = new Int_t*[TrackerParameter.TrackingOrderLines];
        
        Char_t field[100];
        Int_t  aNbOfPlanes;
        for(Int_t iOrder=1 ; iOrder <=TrackerParameter.TrackingOrderLines ; iOrder++)
        {
          nextField();
          sprintf(field, "TrackingOrderLine%d",iOrder);
          cout << "      -> Reading TrackingOrderLine " << iOrder << endl ;
          if( !strcmp( fFieldName, field ))
          {
            nextItem('[');// pre-track planes #
            read_item(aNbOfPlanes);
            nextItem(']');
            cout << "        -> Planes for pre-track = " << aNbOfPlanes << " : { ";
            TrackerParameter.TrackingOrderPreTrack[iOrder-1] = new Int_t[aNbOfPlanes+1];
            TrackerParameter.TrackingOrderPreTrack[iOrder-1][0] = aNbOfPlanes;
            
            nextItem('{');// pre-track planes enumeration starts
            for (Int_t iOrder2=1 ; iOrder2<=aNbOfPlanes ; iOrder2++)
            {
              read_item(TrackerParameter.TrackingOrderPreTrack[iOrder-1][iOrder2]);
              cout << TrackerParameter.TrackingOrderPreTrack[iOrder-1][iOrder2] << " ";
            }
            nextItem('}');// pre-track planes enumeration finished
            cout << "}" << endl;
            
            nextItem('[');// ful-track planes #
            read_item(aNbOfPlanes);
            nextItem(']');
            cout << "        -> Planes for ext-track = " << aNbOfPlanes << " : { ";
            TrackerParameter.TrackingOrderExtTrack[iOrder-1] = new Int_t[aNbOfPlanes+1];
            TrackerParameter.TrackingOrderExtTrack[iOrder-1][0] = aNbOfPlanes;
            
            nextItem('{');// ful-track planes enumeration starts
            for (Int_t iOrder2=1 ; iOrder2<=aNbOfPlanes ; iOrder2++)
            {
              read_item(TrackerParameter.TrackingOrderExtTrack[iOrder-1][iOrder2]);
              cout << TrackerParameter.TrackingOrderExtTrack[iOrder-1][iOrder2] << " ";
            }
            nextItem('}');// pre-track planes enumeration finished            
            cout << "}" << endl ;           
          }
          else
          {
            cout << "        -> Error : can't find line " << field << " in config file" << endl;
            gApplication->Terminate();
          }
        }
	cout << "      -> Reading finished ! " << endl ;
      }
      else
      {
        cout << "      -> ERROR: The field 'TrackingOrderLines' must be >0 ! " << endl;
        gApplication->Terminate();
      }   
    } // end if tracking order

    
    // -*-*-*- Subtrack -*-*-*-
    // JB 2014/12/23
    else if( ! strcmp( fFieldName, "SubtrackNplanes" ) || ! strcmp( fFieldName, "SubtrackNPlanes" ) ) { // if subtrack
      read_item(TrackerParameter.SubtrackNplanes);  
      if ( TrackerParameter.SubtrackNplanes>0 ) {
        TrackerParameter.SubtrackPlanes = new Int_t[TrackerParameter.SubtrackNplanes];
        cout << "  Creating subtrack with " << TrackerParameter.SubtrackNplanes << " planes requested" << endl << "   -> list of planes for subtrack: ";
        nextField();
        if ( ! strcmp( fFieldName, "SubtrackPlanes" ) ) {
          for ( Int_t iplane=0; iplane<TrackerParameter.SubtrackNplanes; iplane++) {
            if( iplane!=0 ) nextItem(':');
            read_item(TrackerParameter.SubtrackPlanes[iplane]);
            cout << TrackerParameter.SubtrackPlanes[iplane] << " ";
          }
          cout << endl;
        }
        else {
          cout << " WARNING: did not detect 'SubtrackPlanes' field so switching SubtrackNplanes to 0!" << endl;
          TrackerParameter.SubtrackNplanes = 0;
        }
      }
    } // end if subtrack


    // ----------------------------------
    //     Vertexing parameters
    // ----------------------------------
    else if( ! strcmp( fFieldName, "VertexMaximum" ) ) { // JB 2013/06/11
      read_item(TrackerParameter.VertexMaximum);  
    }
    else if( ! strcmp( fFieldName, "VertexConstraint" ) ) { // JB 2013/06/21
      read_item(TrackerParameter.VertexConstraint);  
    }
    // ----------------------------------
    //     Alignement parameters
    // ----------------------------------
    else if( ! strcmp( fFieldName, "EventsForAlignmentFit" ) ) {
      read_item(TrackerParameter.EventsForAlignmentFit);  
    }
    else
    {
      if (strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID" ))
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }    
    nextField();
    
  } while ( strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID" ) &&  strcmp( fFieldName, "GeometryName" ));
  
  //Check the correct setup of the experimental parameters
  CheckExpSetupParameters();
  
  if(TrackerParameter.FillNonSensitiveBranch) cout << " Filling up the Non-sensitive branch of output tree!!!" << endl;
  
  cout << endl;
  cout << "Bfield parameters" << endl;
  if(TrackerParameter.BFieldMagnitude == 0.0) cout << "    BFieldMagnitude    parameter not specified. Setting it to defaul value: " << TrackerParameter.BFieldMagnitude << " Tesla" << endl;
  else {
    if(TrackerParameter.BFieldMagnitude < 0) {
      TrackerParameter.BFieldMagnitude = TMath::Abs(TrackerParameter.BFieldMagnitude);
      cout << "    BFieldMagnitude    parameter set to user value is negative. Using its absolute value: " << TrackerParameter.BFieldMagnitude << " Tesla" << endl;
    }
    else cout << "    BFieldMagnitude    parameter set to user value: " << TrackerParameter.BFieldMagnitude << " Tesla" << endl;
  }
  if(TrackerParameter.BFieldDirection == G4ThreeVector(0.0,0.0,1.0)) {
    cout << "    BFieldDirection    parameter not specified. Setting it to defaul value: (" << TrackerParameter.BFieldDirection(0) << "," << TrackerParameter.BFieldDirection(1) << "," << TrackerParameter.BFieldDirection(2) << ")" << endl;
  }
  else {
    //Normalizing Bfield direction vector to unit module
    TrackerParameter.BFieldDirection /= TrackerParameter.BFieldDirection.mag();
    cout << "    BFieldDirection    parameter set to user value: (" << TrackerParameter.BFieldDirection(0) << "," << TrackerParameter.BFieldDirection(1) << "," << TrackerParameter.BFieldDirection(2) << ")" << endl;
  }
  
  
  cout << endl;
  if(TrackerParameter.MediumMaterial == TString("")) {
    TrackerParameter.MediumMaterial = TString("Vacuum");
    cout << "MediumMaterial     parameter not specified. Setting it to defaul value:  " << TrackerParameter.MediumMaterial.Data() << endl;
  }
  else {
    cout << "MediumMaterial     parameter set to user value: " << TrackerParameter.MediumMaterial.Data() << endl;
  }
  cout << endl;
  
  cout << "    -> Reading finished ! " << endl ;
  
  // *******************************************************
  // Test of MANDATORY fields settings // VR 20014/06/30 
  // *******************************************************
  if (TrackerParameter.HitsInPlaneMaximum == -1)
  {
    cout << " ERROR: The field 'HitsInPlaneMaximum' is MANDATORY in config file (section tracker parameters)" << endl;
    gApplication->Terminate();
  }
  
  if(TrackerParameter.TracksMaximum > 0)
  {
    if (TrackerParameter.TracksFinder == 0)
    {
      if (TrackerParameter.SearchHitDistance == -1)
      {
	cout << " ERROR: The field 'SearchHitDistance' is MANDATORY in config file (section tracker parameters)" << endl;
	gApplication->Terminate();
      }
      if (TrackerParameter.PlanesForTrackMinimum == -1)
      {
        cout << " ERROR: The field 'PlanesForTrackMinimum' is MANDATORY in config file (section tracker parameters)" << endl;
        gApplication->Terminate();
      }  
      if (TrackerParameter.HitsInPlaneTrackMaximum == -1)
      {
        cout << " ERROR: The field 'HitsInPlaneTrackMaximum' is MANDATORY in config file (section tracker parameters)" << endl;
        gApplication->Terminate();
      }  
    }
    else if (TrackerParameter.TracksFinder == 1)
    {
      if (TrackerParameter.PlanesForTrackMinimum == -1)
      {
        cout << " ERROR: The field 'PlanesForTrackMinimum' is MANDATORY in config file (section tracker parameters)" << endl;
        gApplication->Terminate();
      }  
      if (TrackerParameter.HitsInPlaneTrackMaximum == -1)
      {
        cout << " ERROR: The field 'HitsInPlaneTrackMaximum' is MANDATORY in config file (section tracker parameters)" << endl;
        gApplication->Terminate();
      }  
      if (TrackerParameter.SearchHitDistance == -1)
      {
	cout << " ERROR: The field 'SearchHitDistance' is MANDATORY in config file (section tracker parameters)" << endl;
	gApplication->Terminate();
      }      
      if (TrackerParameter.SearchMoreHitDistance == -1)
      {
        cout << " ERROR: The field 'SearchMoreHitDistance' is MANDATORY in config file (section tracker parameters) when using TracksFinder: 1" << endl;
        gApplication->Terminate();
      }
    }
    else if (TrackerParameter.TracksFinder == 2)
    {
      if (TrackerParameter.TrackingPass == -1)
      {
        cout << " ERROR: The field 'TrackingPass' is MANDATORY in config file (section tracker parameters) when using TracksFinder: 2" << endl;
        gApplication->Terminate();
      }
      if (TrackerParameter.TrackingOrderLines == -1)
      {
        cout << " ERROR: The field 'TrackingOrderLines' is MANDATORY in config file (section tracker parameters) when using TracksFinder: 2" << endl;
        gApplication->Terminate();
      }
    }
  }

  if(MimosaSimuSetupDebug) {
    cout << " #planes = " << TrackerParameter.Planes << endl;      
    if ( TrackerParameter.Ladders>0 ) {
      cout << " #ladders = " << TrackerParameter.Ladders << endl;
    }
  }

}
//==================================================================
void MimosaSimuSetup::ReadExperimentGeometryParameters() 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //  Experiment Geometrical Parameters
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // VR 2015/01/16
  
  cout << endl << " - Reading Experiment Geometrical Parameters " << endl;
  
  // ################################################
  //            Set default values 
  // ################################################
  //Float_t initValue[3] = {-666.0, -666.0, -666.0};
  G4ThreeVector initValue(-666.0, -666.0, -666.0);
  
  sprintf(IviGeometryParameter.GeometryName,"not defined");
  sprintf(IviGeometryParameter.GeometryVersion,"-1");
  IviGeometryParameter.BeamOrigin             = initValue;
  IviGeometryParameter.BeamSlope              = initValue;
  IviGeometryParameter.BeamDisplayStrongBegin = initValue;
  IviGeometryParameter.BeamDisplayStrongStop  = initValue;
  IviGeometryParameter.BeamDisplayMediumBegin = initValue;
  IviGeometryParameter.BeamDisplayMediumStop  = initValue;
  IviGeometryParameter.BeamDisplayLightBegin  = initValue;
  IviGeometryParameter.BeamDisplayLightStop   = initValue;
  sprintf(IviGeometryParameter.TargetType,"not defined");
  IviGeometryParameter.TargetSize             = initValue;
  IviGeometryParameter.TargetRadius = -1.;
  IviGeometryParameter.TargetLength = -1.;
  //sprintf(IviGeometryParameter.TargetAxis,"O");
  IviGeometryParameter.TargetAxis = TString("O");
  IviGeometryParameter.TargetCenter           = initValue;
  IviGeometryParameter.TrackerOrigin          = initValue;
  IviGeometryParameter.TrackerTilt            = initValue;
  sprintf(IviGeometryParameter.VertexingMethod,"not defined");

  if ( !strcmp(fFieldName,"Inputs") ||  !strcmp(fFieldName, "LadderID") )
  {
    cout << "    -> No Experiment Geometrical Parameters to read (next item is '" << fFieldName << "')" << endl ;
    return ;
  }

  // ----------------------------------
  //     IVI
  // ----------------------------------
  cout << "    -> IVI geometry : " << endl ;
  do
  {  
    if ( ! strcmp( fFieldName, "GeometryName") ) 
    {
      read_strings(IviGeometryParameter.GeometryName, IviGeometryParameter.tpsz);  
      cout << "      -> GeometryName:          " << IviGeometryParameter.GeometryName << endl ;      
    }
    else if ( ! strcmp( fFieldName, "GeometryVersion") ) 
    {
      read_strings(IviGeometryParameter.GeometryVersion, IviGeometryParameter.tpsz);  
      cout << "      -> GeometryVersion:       " << IviGeometryParameter.GeometryVersion << endl ;      
    }
    else if ( ! strcmp( fFieldName, "BeamOriginI") ) 
    {
      read_r3(IviGeometryParameter.BeamOrigin); 
      cout << "      -> BeamOrigin [mm]:    I: " << IviGeometryParameter.BeamOrigin(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.BeamOrigin(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.BeamOrigin(2) << endl ;
    }
    else if ( ! strcmp( fFieldName, "BeamSlopeI") ) 
    {
      read_r3(IviGeometryParameter.BeamSlope); 
      cout << "      -> BeamSlope [/K]  :   I: " << IviGeometryParameter.BeamSlope(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.BeamSlope(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.BeamSlope(2) << endl ;
    }
    else if ( ! strcmp( fFieldName, "BeamDisplayStrongBeginI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayStrongBegin); 
      cout << "      -> BeamDisplayStrongBegin [mm]: I: " << IviGeometryParameter.BeamDisplayStrongBegin(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayStrongBegin(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayStrongBegin(2) << endl ;
    }  
    else if ( ! strcmp( fFieldName, "BeamDisplayStrongStopI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayStrongStop); 
      cout << "      -> BeamDisplayStrongStop [mm]:  I: " << IviGeometryParameter.BeamDisplayStrongStop(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayStrongStop(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayStrongStop(2) << endl ;
    } 
    else if ( ! strcmp( fFieldName, "BeamDisplayMediumBeginI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayMediumBegin); 
      cout << "      -> BeamDisplayMediumBegin [mm]: I: " << IviGeometryParameter.BeamDisplayMediumBegin(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayMediumBegin(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayMediumBegin(2) << endl ;
    }  
    else if ( ! strcmp( fFieldName, "BeamDisplayMediumStopI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayMediumStop); 
      cout << "      -> BeamDisplayMediumStop [mm]:  I: " << IviGeometryParameter.BeamDisplayMediumStop(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayMediumStop(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayMediumStop(2) << endl ;
    }   
    else if ( ! strcmp( fFieldName, "BeamDisplayLightBeginI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayLightBegin); 
      cout << "      -> BeamDisplayLightBegin [mm]:  I: " << IviGeometryParameter.BeamDisplayLightBegin(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayLightBegin(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayLightBegin(2) << endl ;
    }  
    else if ( ! strcmp( fFieldName, "BeamDisplayLightStopI") ) 
    {
      read_r3(IviGeometryParameter.BeamDisplayLightStop); 
      cout << "      -> BeamDisplayLightStop [mm]:   I: " << IviGeometryParameter.BeamDisplayLightStop(0) << endl ;
      cout << "                                      J: " << IviGeometryParameter.BeamDisplayLightStop(1) << endl ;
      cout << "                                      K: " << IviGeometryParameter.BeamDisplayLightStop(2) << endl ;
    }  
    else if ( ! strcmp( fFieldName, "TargetType") ) 
    {
      read_strings(IviGeometryParameter.TargetType, IviGeometryParameter.tpsz);  
      cout << "      -> TargetType:            " << IviGeometryParameter.TargetType << endl ;      
    }
    else if ( ! strcmp( fFieldName, "TargetSizeI") ) 
    {
      read_r3(IviGeometryParameter.TargetSize); 
      cout << "      -> TargetSize [mm]  :  I: " << IviGeometryParameter.TargetSize(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.TargetSize(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.TargetSize(2) << endl ;
    }
    else if ( ! strcmp( fFieldName, "TargetRadius") ) 
    {
      read_item(IviGeometryParameter.TargetRadius); 
      cout << "      -> TargetRadius [mm] :    " << IviGeometryParameter.TargetRadius << endl ;
    }
    else if ( ! strcmp( fFieldName, "TargetLength") ) 
    {
      read_item(IviGeometryParameter.TargetLength); 
      cout << "      -> TargetLength [mm] :    " << IviGeometryParameter.TargetLength << endl ;
    }
    else if ( ! strcmp( fFieldName, "TargetAxis") ) 
    {
      Char_t mmm[1];
      read_strings(mmm,1); 
      IviGeometryParameter.TargetAxis = TString(mmm);
      cout << "      -> TargetAxis :           " << IviGeometryParameter.TargetAxis.Data() << endl ;
    }
    else if ( ! strcmp( fFieldName, "TargetCenterI") ) 
    {
      read_r3(IviGeometryParameter.TargetCenter); 
      cout << "      -> TargetCenter [mm]:  I: " << IviGeometryParameter.TargetCenter(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.TargetCenter(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.TargetCenter(2) << endl ;
    }
    else if ( ! strcmp( fFieldName, "TrackerOriginI") ) 
    {
      read_r3(IviGeometryParameter.TrackerOrigin); 
      cout << "      -> TrackerOrigin [mm]: I: " << IviGeometryParameter.TrackerOrigin(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.TrackerOrigin(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.TrackerOrigin(2) << endl ;
    }
    else if ( ! strcmp( fFieldName, "TrackerTiltI") ) 
    {
      read_r3(IviGeometryParameter.TrackerTilt); 
      cout << "      -> TrackerTilt [deg]:  I: " << IviGeometryParameter.TrackerTilt(0) << endl ;
      cout << "                             J: " << IviGeometryParameter.TrackerTilt(1) << endl ;
      cout << "                             K: " << IviGeometryParameter.TrackerTilt(2) << endl ;
    }  
    else if ( ! strcmp( fFieldName, "VertexingMethod") ) 
    {
      read_strings(IviGeometryParameter.VertexingMethod, IviGeometryParameter.tpsz);  
      cout << "      -> VertexingMethod:       " << IviGeometryParameter.VertexingMethod << endl ;      
    }
    
    nextField();
  } while ( strcmp(fFieldName,"Inputs") &&  strcmp(fFieldName,"LadderID") );
  
  cout << "    -> Reading finished ! (next item is '" << fFieldName << "')" << endl ;
  
}
//==================================================================
void MimosaSimuSetup::ReadLadderParameters( Int_t aLadderNumber) 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of a Ladder 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  
  G4ThreeVector initValue(-666.0, -666.0, -666.0);
  
  pLadderParameter[aLadderNumber].LadderID   = -1;
  pLadderParameter[aLadderNumber].Status     = -1;
  sprintf(pLadderParameter[aLadderNumber].Name,"not defined");
  pLadderParameter[aLadderNumber].Planes     = 0;
  pLadderParameter[aLadderNumber].Position   = initValue;
  pLadderParameter[aLadderNumber].Tilt       = initValue;
  pLadderParameter[aLadderNumber].PlaneList  = NULL;
  pLadderParameter[aLadderNumber].PlaneShift = NULL;
  pLadderParameter[aLadderNumber].PlaneTilt  = NULL;
  
  if ( aLadderNumber<0 || TrackerParameter.Ladders<=aLadderNumber ) {
    printf( "WARNING in ReadLadderParameters: trying to add ladder number %d, outside the limits [1 - %d]\n  --> nothing done!\n", aLadderNumber+1, TrackerParameter.Ladders);
    return;
  }
  else {
    cout << endl << " - Reading Parameters of the Ladder " << aLadderNumber+1 << endl;
  }
  
  // Init
  Int_t nbOfPlanesAssociated = 0;
  
  do {
    
    if ( ! strcmp( fFieldName, "LadderID" ) ) {
      read_item(pLadderParameter[aLadderNumber].LadderID);
      if (MimosaSimuSetupDebug) cout << " ID = " << pLadderParameter[aLadderNumber].LadderID << endl;
    }
    else if( ! strcmp( fFieldName, "LadderName" ) ) {
      read_strings( pLadderParameter[aLadderNumber].Name, pLadderParameter[aLadderNumber].tpsz);
      if (MimosaSimuSetupDebug) cout << " Name = " << pLadderParameter[aLadderNumber].Name << endl;    
    }
    else if( ! strcmp( fFieldName, "Status" ) ) {
      read_item( pLadderParameter[aLadderNumber].Status);
      /*if (MimosaSimuSetupDebug)*/ cout << " Status = " << pLadderParameter[aLadderNumber].Status << endl;
    }
    else if( ! strcmp( fFieldName, "LadderPositionX" ) ) {
      read_r3(pLadderParameter[aLadderNumber].Position);
      pLadderParameter[aLadderNumber].Position *= 1000.; // move from millimeters to microns
      if(MimosaSimuSetupDebug) cout << "  Position X = " << pLadderParameter[aLadderNumber].Position[0] << ",  Y = " << pLadderParameter[aLadderNumber].Position[1] << ", Z = " << pLadderParameter[aLadderNumber].Position[2] << " microns" << endl;
    }
    else if( ! strcmp( fFieldName, "LadderTiltZ" ) ) {
      read_r3(pLadderParameter[aLadderNumber].Tilt);
      pLadderParameter[aLadderNumber].Tilt *= TMath::Pi()/180.; // conversion degree to radian
      if(MimosaSimuSetupDebug) cout << "  Tilt X = " << pLadderParameter[aLadderNumber].Tilt[0] << ",  Y = " << pLadderParameter[aLadderNumber].Tilt[1] << ", Z = " << pLadderParameter[aLadderNumber].Tilt[2] << " rad" << endl;
    }
    else if ( ! strcmp( fFieldName, "NbOfPlanes" ) ) {
      read_item(pLadderParameter[aLadderNumber].Planes);
      pLadderParameter[aLadderNumber].PlaneList = new Int_t[pLadderParameter[aLadderNumber].Planes];
      pLadderParameter[aLadderNumber].PlaneShift = new G4ThreeVector[pLadderParameter[aLadderNumber].Planes];
      pLadderParameter[aLadderNumber].PlaneTilt = new G4ThreeVector[pLadderParameter[aLadderNumber].Planes];
      if (MimosaSimuSetupDebug) cout << " Nb of planes = " << pLadderParameter[aLadderNumber].Planes << endl;
    }
    else if( ! strcmp( fFieldName, "IncludedPlane" ) ) {
      if ( nbOfPlanesAssociated >= pLadderParameter[aLadderNumber].Planes ) {
        printf( "WARNING in ReadLadderParameters: trying to associate plane while maximum # planes (%d) already reached --> nothing done!\n", pLadderParameter[aLadderNumber].Planes);
      }
      else { // if planes can still be associated
        read_item(pLadderParameter[aLadderNumber].PlaneList[nbOfPlanesAssociated]);
        if(MimosaSimuSetupDebug) cout << "    Associated plane " << nbOfPlanesAssociated << " is " << pLadderParameter[aLadderNumber].PlaneList[nbOfPlanesAssociated] << endl;
        
        // Force the reading of the two next info
        nextField();
        if ( ! strcmp( fFieldName, "PlaneShiftU" ) ) {
          read_r3(pLadderParameter[aLadderNumber].PlaneShift[nbOfPlanesAssociated]);
          pLadderParameter[aLadderNumber].PlaneShift[nbOfPlanesAssociated] *= 1000.; // conversion from mm to um
          if(MimosaSimuSetupDebug) cout << "    Plane shift U = " << pLadderParameter[aLadderNumber].PlaneShift[nbOfPlanesAssociated][0] << ",  V = " << pLadderParameter[aLadderNumber].PlaneShift[nbOfPlanesAssociated][1] << ", W = " << pLadderParameter[aLadderNumber].PlaneShift[nbOfPlanesAssociated][2] << " microns" << endl;
        }
        else {
          printf( "WARNING in ReadLadderParameters: expecting Plane %d Shift information!\n", nbOfPlanesAssociated);
        }
        
        nextField();
        if ( ! strcmp( fFieldName, "PlaneTiltU" ) ) {
          read_r3(pLadderParameter[aLadderNumber].PlaneTilt[nbOfPlanesAssociated]);
          pLadderParameter[aLadderNumber].PlaneTilt[nbOfPlanesAssociated] *= M_PI/180.; // conversion degree to radian
          if(MimosaSimuSetupDebug) cout << "    Plane Tilt U = " << pLadderParameter[aLadderNumber].PlaneTilt[nbOfPlanesAssociated][0] << ",  V = " << pLadderParameter[aLadderNumber].PlaneTilt[nbOfPlanesAssociated][1] << ", W = " << pLadderParameter[aLadderNumber].PlaneTilt[nbOfPlanesAssociated][2] << " [rad]" << endl;
        }
        else {
          printf( "WARNING in ReadLadderParameters: expecting Plane %d Tilt information!\n", nbOfPlanesAssociated);
        }
        
        nbOfPlanesAssociated++;
      }  // end if planes can still be associated
    } // end if field is IncludedPlane
    else
    {
      if (  strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID") &&  strcmp( fFieldName, "FileHeaderSize") )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }    
        
    nextField();
    
  } while (  strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID") &&  strcmp( fFieldName, "FileHeaderSize") );
  
  if( nbOfPlanesAssociated < pLadderParameter[aLadderNumber].Planes ) {
    printf( "WARNING in ReadLadderParameters: only %d planes have been associated to ladder %d while %d were expected!\n", nbOfPlanesAssociated, pLadderParameter[aLadderNumber].LadderID, pLadderParameter[aLadderNumber].Planes);
  }
  
  fAddedLadders++;

  /*
  if(DSetupDebug) {
    printf("Parameters for ladder %d of type %s:\n", pLadderParameter[aLadderNumber].LadderID, pLadderParameter[aLadderNumber].Name);
    printf("  status = %d\n", pLadderParameter[aLadderNumber].Status);
    printf("  contains %d ladders:\n", pLadderParameter[aLadderNumber].Planes);
    for ( Int_t iplane=0; iplane<pLadderParameter[aLadderNumber].Planes; iplane++) {
      printf("    plane %d -> absolute nb = %d\n", iplane, pLadderParameter[aLadderNumber].PlaneList[iplane]);
      printf ("   relative tilt : %.2f %.2f %.2f \n"
              , pLadderParameter[aLadderNumber].PlaneTilt[iplane](0)
              , pLadderParameter[aLadderNumber].PlaneTilt[iplane](1)
              , pLadderParameter[aLadderNumber].PlaneTilt[iplane](2)); 
      printf ("   relative shift : %.2f %.2f %.2f \n"
              , pLadderParameter[aLadderNumber].PlaneShift[iplane](0)
              , pLadderParameter[aLadderNumber].PlaneShift[iplane](1)
              , pLadderParameter[aLadderNumber].PlaneShift[iplane](2)); 
    }
  }
   */
  
  if ( ! strcmp( fFieldName, "LadderID") ) { // new ladder
    ReadLadderParameters( fAddedLadders );
  } //end if new plane
  else if ( ! strcmp( fFieldName, "Inputs") ) { // new plane
    ReadPlaneParameters( fAddedPlanes );
  } //end if new plane
    
}
//==================================================================
void MimosaSimuSetup::ReadPlaneParameters( Int_t aPlaneNumber) 
{
 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Planes 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified JB 2013/07/17 new parameters HitFinder, IfDigitize, ...
  // Modified JB 2013/08/13 set default valus for MimosaType, Max/MinNStrips
  // Modified JB 2013/08/15 new ChannelOffset, PlaneResolution parameters 
  // Modified JB 2013/08/29 new Digitizer parameters
  // Modified JB 2014/01/07 default for InitialPedestal and InitialNoise
  // Modified JB 2014/04/21 new deformation parameters
  // Modified AP 2014/07/31 new variables for hot pixel list: file name and fake rate cut
  // Modified AP 2014/12 ?  new Resolution options (U,V or by histograms) 

  if ( aPlaneNumber<0 || TrackerParameter.Planes<=aPlaneNumber ) {
    printf( "WARNING in ReadPlaneParameters: trying to add plane number %d, outside the limits [1 - %d]\n  --> nothing done!\n", aPlaneNumber+1, TrackerParameter.Planes);
    return;
  }
  else {
    cout << endl << " - Reading Parameters of the Plane " << aPlaneNumber+1 << endl;
  }
  
  G4ThreeVector initValue(-666.0, -666.0, -666.0);
  
  // Initialize some default values for non mandatory parameters
  pPlaneParameter[aPlaneNumber].Inputs       = 0;
  for(int kkk=0;kkk<fMaxModules;kkk++) {
    pPlaneParameter[aPlaneNumber].ModuleType[kkk]    = -1;
    pPlaneParameter[aPlaneNumber].ModuleNumber[kkk]  = -1;
    pPlaneParameter[aPlaneNumber].InputNumber[kkk]   = -1;
    pPlaneParameter[aPlaneNumber].ChannelNumber[kkk] = -1;
    pPlaneParameter[aPlaneNumber].ChannelOffset[kkk] = -1;
    pPlaneParameter[aPlaneNumber].Channels[kkk]      = -1;
  }
  sprintf(pPlaneParameter[aPlaneNumber].Name,"not defined");
  sprintf(pPlaneParameter[aPlaneNumber].Purpose,"not defined");
  pPlaneParameter[aPlaneNumber].Readout      = 0;
  pPlaneParameter[aPlaneNumber].AnalysisMode = 0;
  pPlaneParameter[aPlaneNumber].CacheSize    = 0;
  pPlaneParameter[aPlaneNumber].DistanceZero = 0.0;
  pPlaneParameter[aPlaneNumber].Position = initValue;
  pPlaneParameter[aPlaneNumber].Tilt     = initValue;
  pPlaneParameter[aPlaneNumber].AlignmentU    = -999.0;
  pPlaneParameter[aPlaneNumber].AlignmentV    = -999.0;
  pPlaneParameter[aPlaneNumber].AlignmentTilt = -999.0;
  pPlaneParameter[aPlaneNumber].Size          = initValue;
  pPlaneParameter[aPlaneNumber].Strips        = initValue;
  pPlaneParameter[aPlaneNumber].Pitch         = initValue;
  pPlaneParameter[aPlaneNumber].StripSize.setX(0.0);
  pPlaneParameter[aPlaneNumber].StripSize.setY(0.0);
  pPlaneParameter[aPlaneNumber].StripSize.setZ(0.0);
  pPlaneParameter[aPlaneNumber].UsingTrackerResolution = true;
  pPlaneParameter[aPlaneNumber].Mapping = -1;
  pPlaneParameter[aPlaneNumber].ThreshNeighbourSN = -1;
  pPlaneParameter[aPlaneNumber].ThreshSeedSN = -1;
  pPlaneParameter[aPlaneNumber].MaxNStrips = 1000;
  pPlaneParameter[aPlaneNumber].MinNStrips = 1;
  pPlaneParameter[aPlaneNumber].Status = -1;
  pPlaneParameter[aPlaneNumber].ParentLadderID = -1;
  pPlaneParameter[aPlaneNumber].HitPositionAlgorithm = 1;
  pPlaneParameter[aPlaneNumber].EtaCoefficientsN = 0;
  for(int kkk=0;kkk<20;kkk++) pPlaneParameter[aPlaneNumber].EtaCoefficient[kkk] = -999.0;
  pPlaneParameter[aPlaneNumber].EtaLowLimit  = -999.0;
  pPlaneParameter[aPlaneNumber].EtaHighLimit = -999.0;
  pPlaneParameter[aPlaneNumber].KappaCoefficientsN = 0;
  for(int kkk=0;kkk<20;kkk++) pPlaneParameter[aPlaneNumber].KappaCoefficient[kkk] = -999.0;
  pPlaneParameter[aPlaneNumber].KappaLowLimit  = -999.0;
  pPlaneParameter[aPlaneNumber].KappaHighLimit = -999.0;
  pPlaneParameter[aPlaneNumber].GammaCoefficientsN = 0;
  for(int kkk=0;kkk<20;kkk++) pPlaneParameter[aPlaneNumber].GammaCoefficient[kkk] = -999.0;
  pPlaneParameter[aPlaneNumber].GammaLowLimit  = -999.0;
  pPlaneParameter[aPlaneNumber].GammaHighLimit = -999.0;
  pPlaneParameter[aPlaneNumber].NoisyStripsN = 0;
  for(int kkk=0;kkk<20;kkk++) pPlaneParameter[aPlaneNumber].NoisyStripsIndex[kkk] = -1;
  
  pPlaneParameter[aPlaneNumber].HitFinder = 0;
  pPlaneParameter[aPlaneNumber].MimosaType = 0;
  pPlaneParameter[aPlaneNumber].MaxNStrips = 1000;
  pPlaneParameter[aPlaneNumber].MinNStrips = 1;
  pPlaneParameter[aPlaneNumber].ClusterLimitRadius = -1;
  pPlaneParameter[aPlaneNumber].PlaneResolution  = TrackerParameter.Resolution;
  pPlaneParameter[aPlaneNumber].PlaneResolutionU = -1;
  pPlaneParameter[aPlaneNumber].PlaneResolutionV = -1;
  pPlaneParameter[aPlaneNumber].IfDigitize = 0;
  for(int kkk=0;kkk<fMaxDigitalThresholds;kkk++) pPlaneParameter[aPlaneNumber].DigitizeThresholds[kkk] = -999.0;
  
  pPlaneParameter[aPlaneNumber].InitialPedestal = 0;
  pPlaneParameter[aPlaneNumber].InitialNoise = 0;
  pPlaneParameter[aPlaneNumber].CommonRegions = 1; // JB 2014/12/29
  pPlaneParameter[aPlaneNumber].IfDeformed = 0;
  pPlaneParameter[aPlaneNumber].ClusterLimit(0) = pPlaneParameter[aPlaneNumber].ClusterLimit(1) = pPlaneParameter[aPlaneNumber].ClusterLimit(3) = -1;
  for ( int i=0; i<7; i++) {
    pPlaneParameter[aPlaneNumber].CoeffLegendreU[i] = 0.;
    pPlaneParameter[aPlaneNumber].CoeffLegendreV[i] = 0.;
  }
  pPlaneParameter[aPlaneNumber].PlaneThickness      = -1;            //Added, AP 2015/03/10
  pPlaneParameter[aPlaneNumber].PlaneMaterial       = TString("");   //Added, AP 2015/03/10
  pPlaneParameter[aPlaneNumber].PlaneMetalThickness = 0.0;           //Added, AP 2016/07/07
  pPlaneParameter[aPlaneNumber].PlaneEpiThickness   = 1.0;           //Added, AP 2016/07/07
  pPlaneParameter[aPlaneNumber].TimeLimit           = -1;            //Added, JB 2015/05/26
  
  //Plane digitization parameters
  pPlaneParameter[aPlaneNumber].PlaneDigitization    = TString("");   //Added, AP 2016/07/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeOcc     = -1.0;          //Added, AP 2016/07/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise   = -1.0;          //Added, AP 2016/07/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeCalib   = -1.0;          //Added, AP 2016/07/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeThre    = -1.0;          //Added, AP 2016/07/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeADCbits = -1;            //Added, AP 2016/22/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin  = -999.0;        //Added, AP 2016/22/19
  pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax  = -999.0;        //Added, AP 2016/22/19

  TString  HotPixelMapFile(""); // ROOT file name with fake rate map, Added, AP 2014/07/31
  Float_t  FakeRateCut  = 1.0;  // Single pixel fake rate cut,        Added, AP 2014/07/31
  Region_t TheRegion;
  TheRegion.R_col[0] = -1;
  TheRegion.R_col[1] = -1;
  TheRegion.R_lin[0] = -1;
  TheRegion.R_lin[1] = -1;

  TString ResolutionFile("");

  std::vector<Float_t> _ResolutionList;
  std::vector<Float_t> _ResolutionUList;
  std::vector<Float_t> _ResolutionVList;
  std::vector<Region_t> _RegionList;
  std::vector<TString>  _ResolutionFileList;
  _ResolutionList.clear();
  _ResolutionUList.clear();
  _ResolutionVList.clear();
  _RegionList.clear();
  _ResolutionFileList.clear();
  pPlaneParameter[aPlaneNumber].PlanePerformancesList.clear();
  
  do {
    
    if ( ! strcmp( fFieldName, "Inputs" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].Inputs);
      if (MimosaSimuSetupDebug) cout << "   " << pPlaneParameter[aPlaneNumber].Inputs << " inputs:" << endl;
      // Force the reading in the given order
      for (Int_t tl = 0; tl < pPlaneParameter[aPlaneNumber].Inputs; tl++){
        nextField();
        read_item(pPlaneParameter[aPlaneNumber].ModuleType[tl]);
        nextField();
        read_item(pPlaneParameter[aPlaneNumber].ModuleNumber[tl]);
        nextField();
        read_item(pPlaneParameter[aPlaneNumber].InputNumber[tl]);
        nextField();
        read_item(pPlaneParameter[aPlaneNumber].ChannelNumber[tl]);
        nextField();
        if( ! strcmp( fFieldName, "ChannelOffset" ) ) { // JB 2013/08/15
          read_item(pPlaneParameter[aPlaneNumber].ChannelOffset[tl]);
          nextField();
        }
        else {
          pPlaneParameter[aPlaneNumber].ChannelOffset[tl] = 1;
        }
        read_item(pPlaneParameter[aPlaneNumber].Channels[tl]);
        if (MimosaSimuSetupDebug) cout << "    input " << tl << " of type " << pPlaneParameter[aPlaneNumber].ModuleType[tl] << " number " << pPlaneParameter[aPlaneNumber].ModuleNumber[tl] << " input " << pPlaneParameter[aPlaneNumber].InputNumber[tl] << " channel " << pPlaneParameter[aPlaneNumber].ChannelNumber[tl] << " channels " << pPlaneParameter[aPlaneNumber].Channels[tl] << endl;
      }
    }
    else if( ! strcmp( fFieldName, "StripselUse" ) ) {
      // do nothing but need to exclude it because of later field "Strips"
//      Int_t ALL = 64; 
//      ChannelAllUse = new Int_t[ALL];
//      for (Int_t a = 0; a < ALL; a++)    ChannelAllUse[a] = 0xFFFF;
      getRidOfLine();
    }
    else if( ! strcmp( fFieldName, "Name" ) ) {
      read_strings( pPlaneParameter[aPlaneNumber].Name, pPlaneParameter[aPlaneNumber].tpsz);
    }
    else if( ! strcmp( fFieldName, "Purpose" ) ) {
      read_strings( pPlaneParameter[aPlaneNumber].Purpose, pPlaneParameter[aPlaneNumber].tpsz);
    }
    else if( ! strcmp( fFieldName, "Readout" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].Readout);
      // Set #events for initialization to 0 for zero-suppressed data
      // JB 2014/12/29
      if ( pPlaneParameter[aPlaneNumber].Readout>=100 ) {
        pPlaneParameter[aPlaneNumber].InitialPedestal = 0;
        pPlaneParameter[aPlaneNumber].InitialNoise = 0;
      }
    }
    else if( ! strcmp( fFieldName, "MimosaType" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].MimosaType); //RDM210509
    }
    else if( ! strcmp( fFieldName, "AnalysisMode" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].AnalysisMode);
    }
    else if( ! strcmp( fFieldName, "HitFinder" ) ) { // JB 2013/07/17
      read_item(pPlaneParameter[aPlaneNumber].HitFinder);
    }
    else if( ! strcmp( fFieldName, "HotPixelMapFile" ) ) {  // AP 2014/07/31
      read_TStrings( HotPixelMapFile, 350);
    }
    else if( ! strcmp( fFieldName, "FakeRateCut" ) ) {      // AP 2014/07/31
      read_item(FakeRateCut);
      //read_item(pPlaneParameter[aPlaneNumber].FakeRateCut);
    }
    else if( ! strcmp( fFieldName, "InitialPedestal" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].InitialPedestal);
    }
    else if( ! strcmp( fFieldName, "InitialNoise" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].InitialNoise);
      pPlaneParameter[aPlaneNumber].InitialPedestal = pPlaneParameter[aPlaneNumber].InitialNoise; // added to initialize when only InitialNoise is specified, JB 2014/01/07
    }
    else if( ! strcmp( fFieldName, "CacheSize" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].CacheSize);
      if(MimosaSimuSetupDebug) printf("  CacheSize : %d     InitialPedestal : %d      InitialNoise : %d\n", 
                             pPlaneParameter[aPlaneNumber].CacheSize,
                             pPlaneParameter[aPlaneNumber].InitialPedestal,
                             pPlaneParameter[aPlaneNumber].InitialNoise);
    }
    else if( ! strcmp( fFieldName, "IfDigitize" ) ) { // JB 2013/07/17
      read_item(pPlaneParameter[aPlaneNumber].IfDigitize);
      // Force the reading of the thresholds right after the IfDigitize field, if non-zero
      if( pPlaneParameter[aPlaneNumber].IfDigitize>0 ) { // if digitization required
        if(MimosaSimuSetupDebug) printf( "  Digitization emulation required over %d thresholds:\n", pPlaneParameter[aPlaneNumber].IfDigitize);
        nextField();
        if( ! strcmp( fFieldName, "DigitizeThresholds" ) ) { // check there are the thresholds
          for (Int_t k = 0; k < pPlaneParameter[aPlaneNumber].IfDigitize; k++) {
            if( k>0 ) nextItem(':'); // already positionned for 1st value
            read_item(pPlaneParameter[aPlaneNumber].DigitizeThresholds[k]);
            if(MimosaSimuSetupDebug) printf( "  thre[%d] = %d\n", k, pPlaneParameter[aPlaneNumber].DigitizeThresholds[k]);
          }
        }
        else { // if the thresholds are missing
          Error("MimosaSimuSetup", "In the configuration file, the field 'DigitizeThresholds' is expected right after the 'IfDigitize' one, with %d thresholds !", pPlaneParameter[aPlaneNumber].IfDigitize);
        }
      }  // end if digitization required
    }    
    else if( ! strcmp( fFieldName, "PositionsX" ) ) {
      read_r3(pPlaneParameter[aPlaneNumber].Position);
      pPlaneParameter[aPlaneNumber].Position *= 1000; // move from millimeters to microns
      if(MimosaSimuSetupDebug) cout << "  Position X = " << pPlaneParameter[aPlaneNumber].Position(0) << ",  Y = " << pPlaneParameter[aPlaneNumber].Position(1) << ", Z = " << pPlaneParameter[aPlaneNumber].Position(2) << " microns" << endl;
    }
    else if( ! strcmp( fFieldName, "TiltZ" ) || ! strcmp( fFieldName, "Tilt1" )  ) {
      // test with Tilt1 added to comply with MAF files, JB 2013/08/19
      read_r3(pPlaneParameter[aPlaneNumber].Tilt);
      pPlaneParameter[aPlaneNumber].Tilt *= M_PI/180.; // conversion degree to radian

      //cout << "Tilt : (" 
      //   << pPlaneParameter[aPlaneNumber].Tilt[0]*180.0/M_PI << "," 
      //   << pPlaneParameter[aPlaneNumber].Tilt[1]*180.0/M_PI << "," 
      //   << pPlaneParameter[aPlaneNumber].Tilt[2]*180.0/M_PI << ") deg" 
      //   << endl;
    }
    else if( ! strcmp( fFieldName, "AlignementU" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].AlignmentU);
    }
    else if( ! strcmp( fFieldName, "AlignementV" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].AlignmentV);
    }
    else if( ! strcmp( fFieldName, "AlignementTilt" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].AlignmentTilt);
      pPlaneParameter[aPlaneNumber].AlignmentTilt *= TMath::DegToRad();
      if(MimosaSimuSetupDebug) cout << "  AlignementU " << pPlaneParameter[aPlaneNumber].AlignmentU << ",  AlignementV " << pPlaneParameter[aPlaneNumber].AlignmentV << ",  AlignementTilt [deg] " << pPlaneParameter[aPlaneNumber].AlignmentTilt << endl;
    }
    else if( ! strcmp( fFieldName, "SizeU" ) ) {
      read_r3(pPlaneParameter[aPlaneNumber].Size);
      pPlaneParameter[aPlaneNumber].Size *= 1000;
      if(MimosaSimuSetupDebug) printf ("  Size : %f %f %f \n",pPlaneParameter[aPlaneNumber].Size[0],pPlaneParameter[aPlaneNumber].Size[1],pPlaneParameter[aPlaneNumber].Size[2]); 
      printf( "WARNING: you specify a Size for your strips, which is now ignored.\n The sensitive area is now computed from the Pitch and the nb of strips.\n");
    }
    else if( ! strcmp( fFieldName, "StripsU" ) ) {
      read_r3(pPlaneParameter[aPlaneNumber].Strips);
      if(MimosaSimuSetupDebug) printf ("  #strips : %f %f %f \n",pPlaneParameter[aPlaneNumber].Strips[0],pPlaneParameter[aPlaneNumber].Strips[1],pPlaneParameter[aPlaneNumber].Strips[2]); 
    }
    else if( ! strcmp( fFieldName, "PitchU" ) ) {
      read_r3(pPlaneParameter[aPlaneNumber].Pitch);
      pPlaneParameter[aPlaneNumber].Pitch *= 1000;
      if(MimosaSimuSetupDebug) printf ("  strip pitch : %f %f %f \n",pPlaneParameter[aPlaneNumber].Pitch[0],pPlaneParameter[aPlaneNumber].Pitch[1],pPlaneParameter[aPlaneNumber].Pitch[2]); 
      // If Strip Size is not defined yet, set to Pitch, JB 2014/04/21
      if( pPlaneParameter[aPlaneNumber].StripSize[0]<1. ) {
        pPlaneParameter[aPlaneNumber].StripSize = pPlaneParameter[aPlaneNumber].Pitch;
      }
    }
    else if( ! strcmp( fFieldName, "StripSize" ) || ! strcmp( fFieldName, "StripSizeU" )) {
      read_r3(pPlaneParameter[aPlaneNumber].StripSize);
      pPlaneParameter[aPlaneNumber].StripSize *= 1000;
      if(MimosaSimuSetupDebug) printf ("  strip size : %f %f %f \n",pPlaneParameter[aPlaneNumber].StripSize[0],pPlaneParameter[aPlaneNumber].StripSize[1],pPlaneParameter[aPlaneNumber].StripSize[2]); 
    }
    else if( ! strcmp( fFieldName, "Mapping" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].Mapping);
      if(MimosaSimuSetupDebug) cout << "  Mapping = " << pPlaneParameter[aPlaneNumber].Mapping << endl;      
    }
    else if( ! strcmp( fFieldName, "ThreshNeighbourSN" ) || ! strcmp( fFieldName, "ThresNeighbourSN" ) ) { // additional condition (typo), JB 2013/08/19 
      read_item(pPlaneParameter[aPlaneNumber].ThreshNeighbourSN);
    }
    else if( ! strcmp( fFieldName, "ThreshSeedSN" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].ThreshSeedSN);
      if(MimosaSimuSetupDebug) cout << "  Threshold SNR cuts, neighbour = " << pPlaneParameter[aPlaneNumber].ThreshNeighbourSN << " seed = " << pPlaneParameter[aPlaneNumber].ThreshSeedSN << endl ;
    }
    else if( ! strcmp( fFieldName, "MaxNStrips" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].MaxNStrips);
    }
    else if( ! strcmp( fFieldName, "MinNStrips" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].MinNStrips);
    }
    else if( ! strcmp( fFieldName, "ClusterLimit" ) || ! strcmp( fFieldName, "ClusterLimitU" )) {
      read_r3(pPlaneParameter[aPlaneNumber].ClusterLimit);
      pPlaneParameter[aPlaneNumber].ClusterLimit *= 1000;
    }
    else if( ! strcmp( fFieldName, "ClusterLimitRadius" ) ) { // VR 2014/07/16
      read_item(pPlaneParameter[aPlaneNumber].ClusterLimitRadius);
      pPlaneParameter[aPlaneNumber].ClusterLimitRadius *= 1000;
    }
    else if( ! strcmp( fFieldName, "CommonRegions" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].CommonRegions);
    }
    else if( ! strcmp( fFieldName, "Status" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].Status);
      if(MimosaSimuSetupDebug) cout << "  Status of plane = " << pPlaneParameter[aPlaneNumber].Status << endl;
    }
    else if( ! strcmp( fFieldName, "ParentLadder" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].ParentLadderID);
      if(MimosaSimuSetupDebug) cout << "  Parent ladder ID = " << pPlaneParameter[aPlaneNumber].ParentLadderID << endl;
    }
    else if( ! strcmp( fFieldName, "PositionAlgorithm" ) ) {
      read_item(pPlaneParameter[aPlaneNumber].HitPositionAlgorithm);
    }
    else if( ! strcmp( fFieldName, "PlaneResolution" ) ) { // JB 2013/08/15
      Float_t resolution;
      read_item(resolution);
      _ResolutionList.push_back(resolution);
      //read_item(pPlaneParameter[aPlaneNumber].PlaneResolution);
    }
    else if( ! strcmp( fFieldName, "PlaneResolutionU" ) ) { // AP 2014/11/20
      Float_t resolution;
      read_item(resolution);
      _ResolutionUList.push_back(resolution);
      //read_item(pPlaneParameter[aPlaneNumber].PlaneResolutionU);
    }
    else if( ! strcmp( fFieldName, "PlaneResolutionV" ) ) { // AP 2014/11/20
      Float_t resolution;
      read_item(resolution);
      _ResolutionVList.push_back(resolution);
      //read_item(pPlaneParameter[aPlaneNumber].PlaneResolutionV);
    }
    else if( ! strcmp( fFieldName, "PlaneThickness" ) ) { // AP 2015/03/10
      read_item(pPlaneParameter[aPlaneNumber].PlaneThickness);
    }
    else if( ! strcmp( fFieldName, "PlaneMaterial" ) ) {  // AP 2015/03/10
      read_TStrings( pPlaneParameter[aPlaneNumber].PlaneMaterial, 350);
    }
    else if( ! strcmp( fFieldName, "PlaneMetalThickness" ) ) { // AP 2016/07/07
      read_item(pPlaneParameter[aPlaneNumber].PlaneMetalThickness);
    }
    else if( ! strcmp( fFieldName, "PlaneEpiThickness" ) ) {   // AP 2016/07/07
      read_item(pPlaneParameter[aPlaneNumber].PlaneEpiThickness);
    }
    //Plane digitization parameters
    else if( ! strcmp( fFieldName, "PlaneDigitization" ) ) {  // AP 2016/07/19
      read_TStrings( pPlaneParameter[aPlaneNumber].PlaneDigitization, 350);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeOcc" ) ) {   // AP 2016/07/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeOcc);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeNoise" ) ) { // AP 2016/07/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeCalib" ) ) { // AP 2016/07/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeCalib);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeThre" ) ) {  // AP 2016/07/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeThre);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeADCbits" ) ) {  // AP 2016/22/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeADCbits);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeADCMin" ) ) {  // AP 2016/22/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin);
    }
    else if( ! strcmp( fFieldName, "PlaneDigitizeADCMax" ) ) {  // AP 2016/22/19
      read_item(pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax);
    }
    else if( ! strcmp( fFieldName, "ResolutionFile" ) ) {  // AP 2014/07/31
      read_TStrings( ResolutionFile, 350);
      _ResolutionFileList.push_back(ResolutionFile); 
    }
    else if( ! strcmp( fFieldName, "ResolutionRegion" ) ) { // AP 2014/11/20
      read_item(TheRegion.R_col[0]);
      read_item(TheRegion.R_col[1]);
      read_item(TheRegion.R_lin[0]);
      read_item(TheRegion.R_lin[1]);

      //cout << TheRegion.R_col[0] << "   "
      //   << TheRegion.R_col[1] << "   "
      //   << TheRegion.R_lin[0] << "   "
      //   << TheRegion.R_lin[1] << "   "
      //   << endl;

      if(TheRegion.R_col[0] == -1 || 
	 TheRegion.R_col[1] == -1 || 
	 TheRegion.R_lin[0] == -1 || 
	 TheRegion.R_lin[1] == -1) {
	cout << endl;
	cout << "When specifying a region need to give 4 parameters. col_min, col_max, line_min and line_max. Check your inputs. Exiting now!!!" << endl;
	cout << endl;
	assert(false);
      }

      _RegionList.push_back(TheRegion);

    }
    else if( ! strcmp( fFieldName, "Digitizer" ) || ! strcmp( fFieldName, "Digitize" ) ) { // Digitization parameters, JB 2013/08/29
      read_item(pPlaneParameter[aPlaneNumber].IfDigitize);
      if( pPlaneParameter[aPlaneNumber].IfDigitize>fMaxDigitalThresholds ) {
        Warning("ReadPlaneParameters","Nb of digits specified %d, exceeds limit %d --> set back to limit", pPlaneParameter[aPlaneNumber].IfDigitize, fMaxDigitalThresholds);
        pPlaneParameter[aPlaneNumber].IfDigitize = fMaxDigitalThresholds;
      }
      if(MimosaSimuSetupDebug) cout << "  Digitize over " << pPlaneParameter[aPlaneNumber].IfDigitize << " thresholds:" << endl;
      for ( Int_t ithre=0; ithre<pPlaneParameter[aPlaneNumber].IfDigitize; ithre++) {
        nextItem(':');
        read_item( pPlaneParameter[aPlaneNumber].DigitizeThresholds[ithre]);
        if(MimosaSimuSetupDebug) cout << " " << pPlaneParameter[aPlaneNumber].DigitizeThresholds[ithre] << endl;
      }
      if(MimosaSimuSetupDebug) cout << endl;
    } // end Digitization parameters

    else if( ! strcmp( fFieldName, "Deformed" ) || ! strcmp( fFieldName, "deformed" ) ) { // JB 2014/04/21
      read_item(pPlaneParameter[aPlaneNumber].IfDeformed);
      if(MimosaSimuSetupDebug && pPlaneParameter[aPlaneNumber].IfDeformed>0) cout << "  Deforming plane surface according to parameters:" << endl;
    }
    else if( ! strcmp( fFieldName, "LegendreU" ) || ! strcmp( fFieldName, "legendreU" ) ) { // JB 2014/04/21
      Int_t i=0;
      if(MimosaSimuSetupDebug) cout << "     along U:" ;
      while (i<7) {
        read_item( pPlaneParameter[aPlaneNumber].CoeffLegendreU[i]);
        if(MimosaSimuSetupDebug) cout << " " << pPlaneParameter[aPlaneNumber].CoeffLegendreU[i];
        nextItem(':');
        i++;
      }
      if(MimosaSimuSetupDebug) cout << endl;
    } // end U Deformation parameters
    else if( ! strcmp( fFieldName, "LegendreV" ) || ! strcmp( fFieldName, "legendreV" ) ) { // JB 2014/04/21
      Int_t i=0;
      if(MimosaSimuSetupDebug) cout << "     along V:" ;
      while (i<7) {
        read_item( pPlaneParameter[aPlaneNumber].CoeffLegendreV[i]);
        if(MimosaSimuSetupDebug) cout << " " << pPlaneParameter[aPlaneNumber].CoeffLegendreV[i] << endl;
        nextItem(':');
        i++;
      }
      if(MimosaSimuSetupDebug) cout << endl;
    } // end V Deformation parameters
    else if( ! strcmp( fFieldName, "PlaneTimeLimit" ) || ! strcmp( fFieldName, "planetimelimit" ) ) { // JB 2015/05/26
      read_item(pPlaneParameter[aPlaneNumber].TimeLimit);
    }
    else
    {
      if ( strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID") &&  strcmp( fFieldName, "FileHeaderSize") && strcmp( fFieldName, "ModuleTypes" ) && strcmp( fFieldName, "AcqModuleTypes" ) )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }    

    nextField();
    
  } while ( strcmp( fFieldName, "Inputs") &&  strcmp( fFieldName, "LadderID") &&  strcmp( fFieldName, "FileHeaderSize") && strcmp( fFieldName, "ModuleTypes" ) && strcmp( fFieldName, "AcqModuleTypes" ) );

  
  // Set specific parameters for HitFinder=2
  // AP ?
  
  if (pPlaneParameter[aPlaneNumber].HitFinder == 2)
  {
    if(pPlaneParameter[aPlaneNumber].ClusterLimitRadius == -1)
    {
      cout << " ERROR: The field 'ClusterLimitRadius' is MANDATORY in config file when using HitFinder = 2 (section planes parameters)" << endl;
      gApplication->Terminate();
    }
  }
  if (pPlaneParameter[aPlaneNumber].HitFinder != 2)
  {
    if(pPlaneParameter[aPlaneNumber].ClusterLimit(0) == -1 || pPlaneParameter[aPlaneNumber].ClusterLimit(1) == -1)
    {
      cout << " ERROR: The field 'ClusterLimit' is MANDATORY in config file when using HitFinder != 2 (section planes parameters)" << endl;
      gApplication->Terminate();
    }
  }
  

  //Set of specific parameters for the hit spatial resolution
  //AP 12/01/2015
  if(_ResolutionFileList.size() > 0) {
    pPlaneParameter[aPlaneNumber].UsingTrackerResolution = false;
    //NOTE: currently using histogram with residues (huCGwidth_vs_Mult and hvCGwidth_vs_Mult) instead of resolution. Will correction this later
    TString TheMultName         = TString("npix_c");
    TString TheHistoName_ResolU = TString("huCGwidth_vs_Mult");
    TString TheHistoName_ResolV = TString("hvCGwidth_vs_Mult");
    if(_ResolutionFileList.size() == 1) {
      for(int ireg=0;ireg<int(_ResolutionFileList.size());ireg++) {
        TH1F* hMultiplicity = NULL;
        TH1F* hResolutionU  = NULL;
        TH1F* hResolutionV  = NULL;
        TFile *ResolutionVsMultFile = NULL;
        ResolutionVsMultFile = new TFile(_ResolutionFileList[ireg].Data(),"UPDATE");
        if(ResolutionVsMultFile != NULL) {
          TString TheHistoName;
          
          TheHistoName  = TheMultName;
          hMultiplicity = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hMultiplicity == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          hMultiplicity->Scale(1.0/hMultiplicity->GetEntries());
          
          TheHistoName = TheHistoName_ResolU;
          hResolutionU = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hResolutionU == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          
          TheHistoName = TheHistoName_ResolV;
          hResolutionV = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hResolutionV == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          
          PlanePerformances_t AResolutionVsMult;
          AResolutionVsMult.FakeRate        = 0.0;
          AResolutionVsMult.Region.R_col[0] = 0;
          AResolutionVsMult.Region.R_col[1] = pPlaneParameter[aPlaneNumber].Strips(0);
          AResolutionVsMult.Region.R_lin[0] = 0;
          AResolutionVsMult.Region.R_lin[1] = pPlaneParameter[aPlaneNumber].Strips(1);
          
          AResolutionVsMult.GlobalPlaneResolution  = -1.0;
          AResolutionVsMult.GlobalPlaneResolutionU = -1.0;
          AResolutionVsMult.GlobalPlaneResolutionV = -1.0;
          
          AResolutionVsMult.MultProb.clear();
          AResolutionVsMult.MultProbCumul.clear();
          AResolutionVsMult.ResolutionU.clear();
          AResolutionVsMult.ResolutionV.clear();
          for(int imult=0;imult<hResolutionU->GetXaxis()->GetNbins();imult++) {
            AResolutionVsMult.ResolutionU.push_back(hResolutionU->GetBinContent(imult+1));
            AResolutionVsMult.ResolutionV.push_back(hResolutionV->GetBinContent(imult+1));
          }
          for(int imult=0;imult<hResolutionU->GetXaxis()->GetNbins();imult++) {
            if(imult+1 < hResolutionU->GetXaxis()->GetNbins()) {
              double prob = hMultiplicity->GetBinContent(hMultiplicity->FindBin(imult+1));
              AResolutionVsMult.MultProb.push_back(prob);
            }
            else {
              double prob = 0.0;
              for(int kkk=0;kkk<hMultiplicity->GetXaxis()->GetNbins();kkk++) {
                double c = hMultiplicity->GetBinCenter(kkk+1);
                if(c >= imult+1) {
                  prob += hMultiplicity->GetBinContent(kkk+1);
                }
              }
              AResolutionVsMult.MultProb.push_back(prob);
            }
          }
          double prob_cumul = 0.0;
          for(int imult=0;imult<int(AResolutionVsMult.MultProb.size());imult++) {
            prob_cumul += AResolutionVsMult.MultProb[imult];
            AResolutionVsMult.MultProbCumul.push_back(prob_cumul);
          }
          pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
          
          ResolutionVsMultFile->Close();
        } // end if of ROOT file exists
        else {
          cout << endl;
          cout << "Resolution vs multiplicity ROOT file " << _ResolutionFileList[ireg].Data() 
          << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
          cout << endl;
          assert(false);
        }
      }
    }
    else {
      if(_ResolutionFileList.size() != _RegionList.size()) {
        cout << endl;
        cout << "Number of Resolution files given (" << _ResolutionFileList.size()
        << ") is different from number of Regions given (" << _RegionList.size()
        << "). Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      //Check that the regions defined agree with the sensor limits and that they not overlap with each others
      bool RegionsOverlap = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        for(int jreg=ireg+1;jreg<int(_RegionList.size());jreg++) {
          int Corners[4][2];
          Corners[0][0] = _RegionList[jreg].R_col[0]; Corners[0][1] = _RegionList[jreg].R_lin[0];
          Corners[1][0] = _RegionList[jreg].R_col[1]; Corners[1][1] = _RegionList[jreg].R_lin[0];
          Corners[2][0] = _RegionList[jreg].R_col[0]; Corners[2][1] = _RegionList[jreg].R_lin[1];
          Corners[3][0] = _RegionList[jreg].R_col[1]; Corners[3][1] = _RegionList[jreg].R_lin[1];
          
          for(int icorners=0;icorners<4;icorners++) {
            if((Corners[icorners][0] >= _RegionList[ireg].R_col[0] && Corners[icorners][0] <= _RegionList[ireg].R_col[1]) && 
               (Corners[icorners][1] >= _RegionList[ireg].R_lin[0] && Corners[icorners][1] <= _RegionList[ireg].R_lin[1])) {
              RegionsOverlap = true;
              cout << "region " << ireg+1 << " overlaps with region " << jreg+1 << endl;
              break;
            }
          }
          if(RegionsOverlap) break;
        }
        if(RegionsOverlap) break;
      }
      if(RegionsOverlap) {
        cout << endl;
        cout << "At least a couple of the regions defined for plane " << aPlaneNumber << " overlap. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      bool BadRegion = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        if(_RegionList[ireg].R_col[0] < 0 || 
           _RegionList[ireg].R_col[1] > pPlaneParameter[aPlaneNumber].Strips(0)-1 || 
           _RegionList[ireg].R_lin[0] < 0 ||
           _RegionList[ireg].R_lin[1] > pPlaneParameter[aPlaneNumber].Strips(1)-1
           ) {
          BadRegion = true;
          cout << "Region " << ireg+1 << " has limits outside sensor limits!. Region col = (" << _RegionList[ireg].R_col[0] << "," << _RegionList[ireg].R_col[1] << "), lin = ("
          << _RegionList[ireg].R_lin[0] << "," << _RegionList[ireg].R_lin[1] << "). Sensor col = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(0)-1 
          << "), lin = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(1)-1 << ")"
          << endl;
          break;
        }
      }
      if(BadRegion) {
        cout << endl;
        cout << "At least a of the regions defined for Plane " << aPlaneNumber+1 << " overlap is outside of the sensor limits. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        for(int jreg=ireg+1;jreg<int(_RegionList.size());jreg++) {
          if(_ResolutionFileList[ireg] == _ResolutionFileList[jreg]) {
            cout << "WARNING: files for regions " << ireg+1 << " and " << jreg+1 << " for Plane " << aPlaneNumber+1 << " are the same!!!." << endl;
          }
        }
      }
      
      for(int ireg=0;ireg<int(_ResolutionFileList.size());ireg++) {
        TH1F* hMultiplicity = NULL;
        TH1F* hResolutionU  = NULL;
        TH1F* hResolutionV  = NULL;
        TFile *ResolutionVsMultFile = NULL;
        ResolutionVsMultFile = new TFile(_ResolutionFileList[ireg].Data(),"UPDATE");
        if(ResolutionVsMultFile != NULL) {
          TString TheHistoName;
          
          TheHistoName  = TheMultName;
          hMultiplicity = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hMultiplicity == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          hMultiplicity->Scale(1.0/hMultiplicity->GetEntries());
          
          TheHistoName = TheHistoName_ResolU;
          hResolutionU = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hResolutionU == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          
          TheHistoName = TheHistoName_ResolV;
          hResolutionV = (TH1F*)ResolutionVsMultFile->Get(TheHistoName.Data());
          if(hResolutionV == NULL) {
            cout << endl;
            cout << "Histogram with name " << TheHistoName.Data() << " is not found in ROOT file " << _ResolutionFileList[ireg].Data() 
            << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
          
          PlanePerformances_t AResolutionVsMult;
          AResolutionVsMult.FakeRate        = 0.0;
          AResolutionVsMult.Region.R_col[0] = _RegionList[ireg].R_col[0];
          AResolutionVsMult.Region.R_col[1] = _RegionList[ireg].R_col[1];
          AResolutionVsMult.Region.R_lin[0] = _RegionList[ireg].R_lin[0];
          AResolutionVsMult.Region.R_lin[1] = _RegionList[ireg].R_lin[1];
          
          AResolutionVsMult.GlobalPlaneResolution  = -1.0;
          AResolutionVsMult.GlobalPlaneResolutionU = -1.0;
          AResolutionVsMult.GlobalPlaneResolutionV = -1.0;
          
          AResolutionVsMult.MultProb.clear();
          AResolutionVsMult.MultProbCumul.clear();
          AResolutionVsMult.ResolutionU.clear();
          AResolutionVsMult.ResolutionV.clear();
          for(int imult=0;imult<hResolutionU->GetXaxis()->GetNbins();imult++) {
            AResolutionVsMult.ResolutionU.push_back(hResolutionU->GetBinContent(imult+1));
            AResolutionVsMult.ResolutionV.push_back(hResolutionV->GetBinContent(imult+1));
          }
          for(int imult=0;imult<hResolutionU->GetXaxis()->GetNbins();imult++) {
            if(imult+1 < hResolutionU->GetXaxis()->GetNbins()) {
              double prob = hMultiplicity->GetBinContent(hMultiplicity->FindBin(imult+1));
              AResolutionVsMult.MultProb.push_back(prob);
            }
            else {
              double prob = 0.0;
              for(int kkk=0;kkk<hMultiplicity->GetXaxis()->GetNbins();kkk++) {
                double c = hMultiplicity->GetBinCenter(kkk+1);
                if(c >= imult+1) {
                  prob += hMultiplicity->GetBinContent(kkk+1);
                }
              }
              AResolutionVsMult.MultProb.push_back(prob);
            }
          }
          double prob_cumul = 0.0;
          for(int imult=0;imult<int(AResolutionVsMult.MultProb.size());imult++) {
            prob_cumul += AResolutionVsMult.MultProb[imult];
            AResolutionVsMult.MultProbCumul.push_back(prob_cumul);
          }
          pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
          
          ResolutionVsMultFile->Close();
        } // end if of ROOT file exists
        else {
          cout << endl;
          cout << "Resolution vs multiplicity ROOT file " << _ResolutionFileList[ireg].Data() 
          << " for plane " << aPlaneNumber << " and region " << ireg+1 << "doesn't exits. Check your inputs. Exiting now!!!" << endl;
          cout << endl;
          assert(false);
        }
      }
      
    }
  }
  else if(_ResolutionUList.size() > 0 || _ResolutionVList.size() > 0) {
    pPlaneParameter[aPlaneNumber].UsingTrackerResolution = false;
    if(_ResolutionUList.size() != _ResolutionVList.size()) {
      cout << endl;
      cout << "Need to specify both the same number of PlaneResolutionU and PlaneResolutionV parameters for Plane " << aPlaneNumber+1 
      << ". Check your inputs. Exiting now!!!"
      << endl;
      cout << endl;
      assert(false);
    }
    
    if(_ResolutionUList.size() == 1 && _ResolutionVList.size() == 1) {
      PlanePerformances_t AResolutionVsMult;
      AResolutionVsMult.FakeRate        = 0.0;
      AResolutionVsMult.Region.R_col[0] = 0;
      AResolutionVsMult.Region.R_col[1] = pPlaneParameter[aPlaneNumber].Strips(0);
      AResolutionVsMult.Region.R_lin[0] = 0;
      AResolutionVsMult.Region.R_lin[1] = pPlaneParameter[aPlaneNumber].Strips(1);
      
      AResolutionVsMult.GlobalPlaneResolution  = -1.0;
      AResolutionVsMult.GlobalPlaneResolutionU = _ResolutionUList[0];
      AResolutionVsMult.GlobalPlaneResolutionV = _ResolutionVList[0];
      
      AResolutionVsMult.MultProb.clear();
      AResolutionVsMult.MultProbCumul.clear();
      AResolutionVsMult.ResolutionU.clear();
      AResolutionVsMult.ResolutionV.clear();
      pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
    }
    else {
      if(_ResolutionUList.size() != _RegionList.size()) {
        cout << endl;
        cout << "Number of PlaneResolutionU/PlaneResolutionV parameters given (" << _ResolutionUList.size()
        << ") is different from number of Regions given (" << _RegionList.size()
        << "). Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      
      //Check that the regions defined agree with the sensor limits and that they not overlap with each others
      bool RegionsOverlap = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        for(int jreg=ireg+1;jreg<int(_RegionList.size());jreg++) {
          int Corners[4][2];
          Corners[0][0] = _RegionList[jreg].R_col[0]; Corners[0][1] = _RegionList[jreg].R_lin[0];
          Corners[1][0] = _RegionList[jreg].R_col[1]; Corners[1][1] = _RegionList[jreg].R_lin[0];
          Corners[2][0] = _RegionList[jreg].R_col[0]; Corners[2][1] = _RegionList[jreg].R_lin[1];
          Corners[3][0] = _RegionList[jreg].R_col[1]; Corners[3][1] = _RegionList[jreg].R_lin[1];
          
          for(int icorners=0;icorners<4;icorners++) {
            if((Corners[icorners][0] >= _RegionList[ireg].R_col[0] && Corners[icorners][0] <= _RegionList[ireg].R_col[1]) && 
               (Corners[icorners][1] >= _RegionList[ireg].R_lin[0] && Corners[icorners][1] <= _RegionList[ireg].R_lin[1])) {
              RegionsOverlap = true;
              cout << "region " << ireg+1 << " overlaps with region " << jreg+1 << endl;
              break;
            }
          }
          if(RegionsOverlap) break;
        }
        if(RegionsOverlap) break;
      }
      
      if(RegionsOverlap) {
        cout << endl;
        cout << "At least a couple of the regions defined for plane " << aPlaneNumber << " overlap. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      bool BadRegion = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        if(_RegionList[ireg].R_col[0] < 0 || 
           _RegionList[ireg].R_col[1] > pPlaneParameter[aPlaneNumber].Strips(0)-1 || 
           _RegionList[ireg].R_lin[0] < 0 ||
           _RegionList[ireg].R_lin[1] > pPlaneParameter[aPlaneNumber].Strips(1)-1
           ) {
          BadRegion = true;
          cout << "Region " << ireg+1 << " has limits outside sensor limits!. Region col = (" << _RegionList[ireg].R_col[0] << "," << _RegionList[ireg].R_col[1] << "), lin = ("
          << _RegionList[ireg].R_lin[0] << "," << _RegionList[ireg].R_lin[1] << "). Sensor col = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(0)-1 
          << "), lin = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(1)-1 << ")"
          << endl;
          break;
        }
      }
      if(BadRegion) {
        cout << endl;
        cout << "At least a of the regions defined for Plane " << aPlaneNumber+1 << " overlap is outside of the sensor limits. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      
      for(int ireg=0;ireg<int(_ResolutionUList.size());ireg++) {
        PlanePerformances_t AResolutionVsMult;
        AResolutionVsMult.FakeRate        = 0.0;
        AResolutionVsMult.Region.R_col[0] = _RegionList[ireg].R_col[0];
        AResolutionVsMult.Region.R_col[1] = _RegionList[ireg].R_col[1];
        AResolutionVsMult.Region.R_lin[0] = _RegionList[ireg].R_lin[0];
        AResolutionVsMult.Region.R_lin[1] = _RegionList[ireg].R_lin[1];
        
        AResolutionVsMult.GlobalPlaneResolution  = -1.0;
        AResolutionVsMult.GlobalPlaneResolutionU = _ResolutionUList[ireg];
        AResolutionVsMult.GlobalPlaneResolutionV = _ResolutionVList[ireg];
        
        AResolutionVsMult.MultProb.clear();
        AResolutionVsMult.MultProbCumul.clear();
        AResolutionVsMult.ResolutionU.clear();
        AResolutionVsMult.ResolutionV.clear();
        pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
      }
      
    }
  }
  else if(_ResolutionList.size() > 0) {
    pPlaneParameter[aPlaneNumber].UsingTrackerResolution = false;
    if(_ResolutionList.size() == 1) {
      PlanePerformances_t AResolutionVsMult;
      AResolutionVsMult.FakeRate        = 0.0;
      AResolutionVsMult.Region.R_col[0] = 0;
      AResolutionVsMult.Region.R_col[1] = pPlaneParameter[aPlaneNumber].Strips(0);
      AResolutionVsMult.Region.R_lin[0] = 0;
      AResolutionVsMult.Region.R_lin[1] = pPlaneParameter[aPlaneNumber].Strips(1);
      
      AResolutionVsMult.GlobalPlaneResolution  = _ResolutionList[0];
      AResolutionVsMult.GlobalPlaneResolutionU = -1;
      AResolutionVsMult.GlobalPlaneResolutionV = -1;
      
      AResolutionVsMult.MultProb.clear();
      AResolutionVsMult.MultProbCumul.clear();
      AResolutionVsMult.ResolutionU.clear();
      AResolutionVsMult.ResolutionV.clear();
      pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
    }
    else {
      if(_ResolutionList.size() != _RegionList.size()) {
        cout << endl;
        cout << "Number of PlaneResolution parameters given (" << _ResolutionUList.size()
        << ") is different from number of Regions given (" << _RegionList.size()
        << "). Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }

      //Check that the regions defined agree with the sensor limits and that they not overlap with each others
      bool RegionsOverlap = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        for(int jreg=ireg+1;jreg<int(_RegionList.size());jreg++) {
          int Corners[4][2];
          Corners[0][0] = _RegionList[jreg].R_col[0]; Corners[0][1] = _RegionList[jreg].R_lin[0];
          Corners[1][0] = _RegionList[jreg].R_col[1]; Corners[1][1] = _RegionList[jreg].R_lin[0];
          Corners[2][0] = _RegionList[jreg].R_col[0]; Corners[2][1] = _RegionList[jreg].R_lin[1];
          Corners[3][0] = _RegionList[jreg].R_col[1]; Corners[3][1] = _RegionList[jreg].R_lin[1];
          
          for(int icorners=0;icorners<4;icorners++) {
            if((Corners[icorners][0] >= _RegionList[ireg].R_col[0] && Corners[icorners][0] <= _RegionList[ireg].R_col[1]) && 
               (Corners[icorners][1] >= _RegionList[ireg].R_lin[0] && Corners[icorners][1] <= _RegionList[ireg].R_lin[1])) {
              RegionsOverlap = true;
              cout << "region " << ireg+1 << " overlaps with region " << jreg+1 << endl;
              break;
            }
          }
          if(RegionsOverlap) break;
        }
        if(RegionsOverlap) break;
      }
      if(RegionsOverlap) {
        cout << endl;
        cout << "At least a couple of the regions defined for Plane " << aPlaneNumber+1 << " overlap. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      bool BadRegion = false;
      for(int ireg=0;ireg<int(_RegionList.size());ireg++) {
        if(_RegionList[ireg].R_col[0] < 0 || 
           _RegionList[ireg].R_col[1] > pPlaneParameter[aPlaneNumber].Strips(0)-1 || 
           _RegionList[ireg].R_lin[0] < 0 ||
           _RegionList[ireg].R_lin[1] > pPlaneParameter[aPlaneNumber].Strips(1)-1
           ) {
          BadRegion = true;
          cout << "Region " << ireg+1 << " has limits outside sensor limits!. Region col = (" << _RegionList[ireg].R_col[0] << "," << _RegionList[ireg].R_col[1] << "), lin = ("
          << _RegionList[ireg].R_lin[0] << "," << _RegionList[ireg].R_lin[1] << "). Sensor col = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(0)-1 
          << "), lin = (" << 0 << "," << pPlaneParameter[aPlaneNumber].Strips(1)-1 << ")"
          << endl;
          break;
        }
      }
      if(BadRegion) {
        cout << endl;
        cout << "At least a of the regions defined for Plane " << aPlaneNumber+1 << " overlap is outside of the sensor limits. Check you inputs. Exiting now!!!!" << endl;
        cout << endl;
        assert(false);
      }
      
      for(int ireg=0;ireg<int(_ResolutionList.size());ireg++) {
        PlanePerformances_t AResolutionVsMult;
        AResolutionVsMult.FakeRate        = 0.0;
        AResolutionVsMult.Region.R_col[0] = _RegionList[ireg].R_col[0];
        AResolutionVsMult.Region.R_col[1] = _RegionList[ireg].R_col[1];
        AResolutionVsMult.Region.R_lin[0] = _RegionList[ireg].R_lin[0];
        AResolutionVsMult.Region.R_lin[1] = _RegionList[ireg].R_lin[1];
        
        AResolutionVsMult.GlobalPlaneResolution  = _ResolutionList[ireg];
        AResolutionVsMult.GlobalPlaneResolutionU = -1.0;
        AResolutionVsMult.GlobalPlaneResolutionV = -1.0;
        
        AResolutionVsMult.MultProb.clear();
        AResolutionVsMult.MultProbCumul.clear();
        AResolutionVsMult.ResolutionU.clear();
        AResolutionVsMult.ResolutionV.clear();
        pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
      }
    }
  }
  else if(_ResolutionList.size() == 0) {
    pPlaneParameter[aPlaneNumber].UsingTrackerResolution = true;
    PlanePerformances_t AResolutionVsMult;
    AResolutionVsMult.FakeRate        = 0.0;
    AResolutionVsMult.Region.R_col[0] = 0;
    AResolutionVsMult.Region.R_col[1] = pPlaneParameter[aPlaneNumber].Strips(0);
    AResolutionVsMult.Region.R_lin[0] = 0;
    AResolutionVsMult.Region.R_lin[1] = pPlaneParameter[aPlaneNumber].Strips(1);
    
    AResolutionVsMult.GlobalPlaneResolution  = TrackerParameter.Resolution;
    AResolutionVsMult.GlobalPlaneResolutionU = -1.0;
    AResolutionVsMult.GlobalPlaneResolutionV = -1.0;
    
    AResolutionVsMult.MultProb.clear();
    AResolutionVsMult.MultProbCumul.clear();
    AResolutionVsMult.ResolutionU.clear();
    AResolutionVsMult.ResolutionV.clear();
    pPlaneParameter[aPlaneNumber].PlanePerformancesList.push_back(AResolutionVsMult);
  }
  
  ///////////////////////////////////
  //Printing out the resolution parameters obtained from config file:
  ///////////////////////////////////
  if(pPlaneParameter[aPlaneNumber].PlanePerformancesList.size() == 1) {
    cout << "Using fake rate of " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].FakeRate << " for the whole matrix" << endl;

    //No region defined:
    if(pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].ResolutionU.size() == 0) {
      //No multiplicity dependent resolution defined
      if(pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].GlobalPlaneResolutionU == -1) {
	//Using same resolution for U and V
	if(pPlaneParameter[aPlaneNumber].UsingTrackerResolution) {
	  //Using same resolution as the one defined on tracker block
	  cout << "Using global TrakerResolution for PlaneResolution for Plane " << aPlaneNumber+1 
	       << ". Resolution = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].GlobalPlaneResolution << "um" << endl;
	  cout << "Using the same resolution for U and V directions!!!" << endl;
	}
	else {
	  //Using user defined global resolution
	  cout << "Using user defined global PlaneResolution for Plane " << aPlaneNumber+1 
	       << ". Resolution = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].GlobalPlaneResolution << "um" << endl;
	  cout << "Using the same resolution for U and V directions!!!" << endl;
	}
      }
      else {
	//Using different global resolutions for U and 
	cout << "Using different PlaneResolutions for Plane " << aPlaneNumber+1 << " for U and V directions." << endl;
	cout << "ResolutionU = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].GlobalPlaneResolutionU << "um" << endl;
	cout << "ResolutionV = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].GlobalPlaneResolutionV << "um" << endl;
      }
    }
    else {
      //Using multiplicity dependent resolutions for U and V
      cout << "Multiplicity dependent resolution for Plane" << aPlaneNumber+1 << ":" << endl;
      for(int imult=0;imult<int(pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].ResolutionU.size());imult++) {
	cout << "(ResolutionU,ResolutionV) = (" 
	     << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].ResolutionU[imult] << ","
	     << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].ResolutionV[imult] << ") um (prob = "
	     << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].MultProb[imult]*100.0 << " %, prob_cumul = "
	     << pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].MultProbCumul[imult]*100.0 << " %) for ";
	if(imult < int(pPlaneParameter[aPlaneNumber].PlanePerformancesList[0].ResolutionU.size()-1)) {
	  cout << " mult = " << imult+1;
	}
	else cout << " mult >= " << imult+1;
	cout << endl;
      }
      cout << endl;
    }
  }
  else {
    //Different regions defined:
    cout << "Found " << pPlaneParameter[aPlaneNumber].PlanePerformancesList.size()
	 << " regions for Plane " << aPlaneNumber+1 << "."
	 << endl;
    for(int ireg=0;ireg<int(pPlaneParameter[aPlaneNumber].PlanePerformancesList.size());ireg++) {
      cout << "Region " << ireg+1 << ": col = (" 
	   << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].Region.R_col[0] << ","
	   << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].Region.R_col[1] << "), lin = (" 
	   << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].Region.R_lin[0] << "," 
	   << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].Region.R_lin[1] << ")."
	   << " Estimated fake rate = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].FakeRate 
	   << endl;
      if(pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].ResolutionU.size() == 0) {
	//No multiplicity dependent resolution defined
	if(pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].GlobalPlaneResolutionU == -1) {
	  //Using same resolution for U and V
	  if(pPlaneParameter[aPlaneNumber].UsingTrackerResolution) {
	    //Using same resolution as the one defined on tracker block
	    cout << "Using global TrakerResolution for PlaneResolution for Plane " << aPlaneNumber+1 << " and Region " << ireg+1
		 << ". Resolution = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].GlobalPlaneResolution << "um" << endl;
	    cout << "Using the same resolution for U and V directions!!!" << endl;
	  }
	  else {
	    //Using user defined global resolution
	    cout << "Using user defined global PlaneResolution for Plane " << aPlaneNumber+1 << " and Region " << ireg+1
		 << ". Resolution = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].GlobalPlaneResolution << "um" << endl;
	    cout << "Using the same resolution for U and V directions!!!" << endl;
	  }
	}
	else {
	  //Using different global resolutions for U and V
	  cout << "Using different PlaneResolutions for Plane " << aPlaneNumber+1 << " for U and V directions." << endl;
	  cout << "ResolutionU = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].GlobalPlaneResolutionU << "um" << endl;
	  cout << "ResolutionV = " << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].GlobalPlaneResolutionV << "um" << endl;
	}
      }
      else {
	//Using multiplicity dependent resolutions for U and V
	for(int imult=0;imult<int(pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].ResolutionU.size());imult++) {
	  cout << "(ResolutionU,ResolutionV) = (" 
	       << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].ResolutionU[imult] << ","
	       << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].ResolutionV[imult] << ") um (prob = "
	       << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].MultProb[imult]*100.0 << " %, prob_cumul = "
	       << pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].MultProbCumul[imult]*100.0 << " %) for ";
	  if(imult < int(pPlaneParameter[aPlaneNumber].PlanePerformancesList[ireg].ResolutionU.size()-1)) {
	    cout << " mult = " << imult+1;
	  }
	  else cout << " mult >= " << imult+1;
	  cout << endl;
	}
	cout << endl;
      }
    }
  }

  if(pPlaneParameter[aPlaneNumber].PlaneThickness < 0) {
    pPlaneParameter[aPlaneNumber].PlaneThickness = 1.0e-10;
    cout << "No specification of PlaneThickness. Use default value of " << pPlaneParameter[aPlaneNumber].PlaneThickness << "um. (Insignificant multiple scattering)" << endl;
  }
  else cout << "PlaneThickness set to user value: " << pPlaneParameter[aPlaneNumber].PlaneThickness << "um." << endl;

  if(pPlaneParameter[aPlaneNumber].PlaneMaterial  == TString("")) {
    pPlaneParameter[aPlaneNumber].PlaneMaterial = TString("silicon");
    cout << "No specification of PlaneMaterial.  Use default value of " << pPlaneParameter[aPlaneNumber].PlaneMaterial.Data() << "." << endl;
  }
  else {
    cout << "PlaneMaterial  set to user value: " << pPlaneParameter[aPlaneNumber].PlaneMaterial.Data() << endl;
  }

  if((pPlaneParameter[aPlaneNumber].PlaneMetalThickness >= 0 && pPlaneParameter[aPlaneNumber].PlaneMetalThickness <= 1) && 
     (pPlaneParameter[aPlaneNumber].PlaneEpiThickness   >= 0 && pPlaneParameter[aPlaneNumber].PlaneEpiThickness   <= 1)
    ) {
    if(pPlaneParameter[aPlaneNumber].PlaneMetalThickness + pPlaneParameter[aPlaneNumber].PlaneEpiThickness > 1) {
      cout << "WARNING : Plane " << aPlaneNumber+1 << " parameters PlaneMetalThickness and PlaneEpiThickness have values which sum exceeds 1. Setting them to default values 0 and 1, respectively!" << endl;
      pPlaneParameter[aPlaneNumber].PlaneMetalThickness = 0.0;
      pPlaneParameter[aPlaneNumber].PlaneEpiThickness   = 1.0;
    }
    else {
      cout << "PlaneMetalThickness set to value: " << pPlaneParameter[aPlaneNumber].PlaneMetalThickness*100 
           << " % of total thickness (" << pPlaneParameter[aPlaneNumber].PlaneMetalThickness*pPlaneParameter[aPlaneNumber].PlaneThickness << " um)" 
	   << endl;
      cout << "PlaneEpiThickness   set to value: " << pPlaneParameter[aPlaneNumber].PlaneEpiThickness*100   
           << " % of total thickness (" << pPlaneParameter[aPlaneNumber].PlaneEpiThickness*pPlaneParameter[aPlaneNumber].PlaneThickness 
           << " um)" << endl;
    }
  }
  else {
    cout << "WARNING : Plane " << aPlaneNumber+1 << " parameters PlaneMetalThickness and/or PlaneEpiThickness have values outside valid range (0,1). Setting them to default values 0 and 1, respectively!" << endl;
    pPlaneParameter[aPlaneNumber].PlaneMetalThickness = 0.0;
    pPlaneParameter[aPlaneNumber].PlaneEpiThickness   = 1.0;
  }

  //Check the correct setup of digitization parameters
  CheckPlaneDigitizationParameters(aPlaneNumber);
  
  //Filling up list of sensitive planes
  if(pPlaneParameter[aPlaneNumber].Readout > 0) ListOfSensitivePlanes.push_back(aPlaneNumber+1);
  
  // ---------
  // This is the end of plane parameters

  fAddedPlanes++;
    
  if ( ! strcmp( fFieldName, "Inputs") ) { // new plane
    ReadPlaneParameters( fAddedPlanes );
  } //end if new plane
  //else if ( ! strcmp( fFieldName, "LadderID") ) { // new ladder
    //ReadLadderParameters( fAddedLadders );
  //} //end if new plane
  
}
//==================================================================
void MimosaSimuSetup::ReadDAQParameters() 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter of the Data Acquisition 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16

  cout << endl << " - Reading  Parameter of the Data Acquisition " << endl;
  
  AcqParameter.FileHeaderSize    = -1;
  AcqParameter.EventBufferSize   = -1;
  AcqParameter.FileHeaderLine    = -1;
  AcqParameter.ModuleTypes       = -1;
  AcqParameter.BinaryCoding      = -1;
  AcqParameter.TriggerMode       = -1;
  AcqParameter.EventBuildingMode =  1; //SS 2011.11.14. Setting EventBuildingMode=1 by default. It can be changed during initialisation of the session or from the config file

  do {
    
    if( ! strcmp( fFieldName, "FileHeaderSize" ) ) {
      read_item(AcqParameter.FileHeaderSize);
    }
    else if( ! strcmp( fFieldName, "EventBufferSize" ) ) {
      read_item(AcqParameter.EventBufferSize);
    }
    else if( ! strcmp( fFieldName, "FileHeaderLine" ) || ! strcmp( fFieldName, "FileHeaderLine[d]" )) {
      read_item(AcqParameter.FileHeaderLine);
    }
    else if( ! strcmp( fFieldName, "ModuleTypes" ) || ! strcmp( fFieldName, "AcqModuleTypes" )  ) {
      read_item(AcqParameter.ModuleTypes);
    }
    else if( ! strcmp( fFieldName, "BinaryCoding" ) ) {
      read_item(AcqParameter.BinaryCoding);
    }
    else if( ! strcmp( fFieldName, "TriggerMode" ) ) {
      read_item(AcqParameter.TriggerMode);
    }
    else if( ! strcmp( fFieldName, "EventBuildingMode" ) ) {
      read_item(AcqParameter.EventBuildingMode);
    }
    else
    {
      if (  strcmp( fFieldName, "Name") )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }      
    nextField();
    
  } while (  strcmp( fFieldName, "Name") );
  
  
  if( MimosaSimuSetupDebug) {
    cout << "   header file size " << AcqParameter.FileHeaderSize << endl;
    cout << "   event buffer size " << AcqParameter.EventBufferSize << endl;
    cout << "   event header size " << AcqParameter.FileHeaderLine << endl;
    cout << "   nb of mod types " << AcqParameter.ModuleTypes << endl;
    cout << "   endian coding " << AcqParameter.BinaryCoding << endl;
    cout << "   trigger mode " << AcqParameter.TriggerMode << endl;
    cout << "   event building mode " << AcqParameter.EventBuildingMode << endl;
  }
  
}
//==================================================================
void MimosaSimuSetup::ReadDAQBoardParameters( Int_t aBoardNumber) 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter of the Data Acquisition boards in this run
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified: JB 2013/06/20 correction to call ReadDAQBoardParameters at the end
  // Modified: JB 2013/06/22 new EventBuildingBoardMode param
  // Modified: JB 2013/08/16 manage different parameters for each input
  // Modified: JB+CB+PLR 2015/03/24 additional parameter for INFN decoder

  
  sprintf(pAcqModuleParameter[aBoardNumber].Name,"");
  pAcqModuleParameter[aBoardNumber].Devices                   = -1;
  pAcqModuleParameter[aBoardNumber].Type                      = -1;
  pAcqModuleParameter[aBoardNumber].EventBuildingBoardMode    = -1;
  pAcqModuleParameter[aBoardNumber].Inputs                    = -1;
  for(int i=0;i<fMaxModules;i++) {
    pAcqModuleParameter[aBoardNumber].Channels[i]             = -1;
    pAcqModuleParameter[aBoardNumber].Bits[i]                 = -1;
    pAcqModuleParameter[aBoardNumber].SigBits[i]              = -1;
    pAcqModuleParameter[aBoardNumber].DeviceDataFile[i]       = NULL;
    pAcqModuleParameter[aBoardNumber].NbOfFramesPerChannel[i] = -1;
    pAcqModuleParameter[aBoardNumber].PixelShiftMod[i]        = -1;
  }
  pAcqModuleParameter[aBoardNumber].FirstTriggerChannel       = -1;
  pAcqModuleParameter[aBoardNumber].LastTriggerChannel        = -1;
  pAcqModuleParameter[aBoardNumber].PixelShift                = -1;
  pAcqModuleParameter[aBoardNumber].AmpOffset                 = -1;
  pAcqModuleParameter[aBoardNumber].AmpFactor                 = -1;
  pAcqModuleParameter[aBoardNumber].Trailer                   = 0;
  
  if ( aBoardNumber<0 || AcqParameter.ModuleTypes<=aBoardNumber ) {
    printf( "WARNING in ReadDAQBoardParameters: trying to add a board type %d, outside the limits [1 - %d]\n  --> nothing done!\n", aBoardNumber+1, AcqParameter.ModuleTypes);
    return;
  }
  else {
    cout << endl << " - Reading Parameters of the Data Acquisition board type " << aBoardNumber+1 << endl;
  }
  
  // Some default values
  pAcqModuleParameter[aBoardNumber].EventBuildingBoardMode = AcqParameter.EventBuildingMode; // default, JB 2013/06/22
  pAcqModuleParameter[aBoardNumber].FirstTriggerChannel = 0;
  pAcqModuleParameter[aBoardNumber].LastTriggerChannel = 0;
  pAcqModuleParameter[aBoardNumber].NbOfFramesPerChannel[0] = 2;
  pAcqModuleParameter[aBoardNumber].PixelShift = 3; // JB,CB,PLR 2015/03/24
  pAcqModuleParameter[aBoardNumber].AmpOffset = 32768; // JB,CB,PLR 2015/03/24
  pAcqModuleParameter[aBoardNumber].AmpFactor = 1.; // JB,CB,PLR 2015/03/24
  pAcqModuleParameter[aBoardNumber].Trailer = 0xfafafafa; // JB,CB,PLR 2015/03/24

  // Taking into account Pixel Shift, which could be defined:
  // - globally with a single PixelShift
  // - for each module with as many PixleShiftMod as module
  for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
    pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod] = -1; // indicate not define
  }
  
  Int_t inputCounter = -1; // starting value for inputs, updated at the next "channels" field
  
  do {
    
    if ( ! strcmp( fFieldName, "Name") ) {
      read_strings( pAcqModuleParameter[aBoardNumber].Name, pAcqModuleParameter[aBoardNumber].tpsz);
    }
    else if( ! strcmp( fFieldName, "Type" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].Type);
      pAcqModuleParameter[aBoardNumber].EventBuildingBoardMode = pAcqModuleParameter[aBoardNumber].Type%10; // JB 2013/06/22
    }
    else if( ! strcmp( fFieldName, "EventBuildingBoardMode" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].EventBuildingBoardMode); // JB 2013/06/22
    }
    else if( ! strcmp( fFieldName, "Devices" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].Devices);
    }
    else if( ! strcmp( fFieldName, "Inputs" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].Inputs);
    }
    else if( ! strcmp( fFieldName, "Channels" ) ) {
      if( inputCounter<=fMaxModules) {
        inputCounter++;  
      }
      else {
        printf("WARNING in DSetup:, trying to add more than %d inputs, ignored\n", fMaxModules);
      }
      read_item(pAcqModuleParameter[aBoardNumber].Channels[inputCounter]);
    }
    else if( ! strcmp( fFieldName, "Bits" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].Bits[inputCounter]);
    }
    else if( ! strcmp( fFieldName, "SignificantBits" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].SigBits[inputCounter]);
    }
    else if( ! strcmp( fFieldName, "FirstTriggerChannel" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].FirstTriggerChannel);
    }
    else if( ! strcmp( fFieldName, "LastTriggerChannel" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].LastTriggerChannel);
    }
    else if( ! strcmp( fFieldName, "NbOfFramesPerChannel" ) ) {
      read_item(pAcqModuleParameter[aBoardNumber].NbOfFramesPerChannel[inputCounter]);
    }
    else if( ! strcmp( fFieldName, "PixelShift" ) && strcmp( fFieldName, "PixelShiftMod" ) ) {  // JB,CB,PLR 2015/03/24
      read_item(pAcqModuleParameter[aBoardNumber].PixelShift);
      // if local PixelShiftMod is not defined, update with global PixelShift
      for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
        if( pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod] == -1 ) {
          pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod] = pAcqModuleParameter[aBoardNumber].PixelShift;
        }
      }

    }
    else if( ! strcmp( fFieldName, "AmpOffset" ) ) { // JB,CB,PLR 2015/03/24
      read_item(pAcqModuleParameter[aBoardNumber].AmpOffset);
    }
    else if( ! strcmp( fFieldName, "AmpFactor" ) ) { // JB,CB,PLR 2015/03/24
      read_item(pAcqModuleParameter[aBoardNumber].AmpFactor);
    }
    else if( ! strcmp( fFieldName, "Trailer" ) ) { // JB,CB,PLR 2015/03/24
      read_item(pAcqModuleParameter[aBoardNumber].Trailer);
      }
    else if( ! strcmp( fFieldName, "DataFile" ) || ! strcmp( fFieldName, "DataFile1" )) {
      // reading data file name for each device of this type, JB 2009/05/25
      for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
        pAcqModuleParameter[aBoardNumber].DeviceDataFile[iMod] = new Char_t[pAcqModuleParameter[aBoardNumber].tpsz];
        pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod] = 3; // default value
        if( iMod>0 ) nextField();
        if( strstr( fFieldName, "DataFile" ) ) {
          read_strings( pAcqModuleParameter[aBoardNumber].DeviceDataFile[iMod], pAcqModuleParameter[aBoardNumber].tpsz);
        }
        else {
          printf( "WARNING in ReadDAQBoardParameters: field %s found while 'DataFile%d' expected!\n    -> Board %d of type %d lacks rawdata connection.", fFieldName, iMod+1, iMod+1, aBoardNumber);
        }
      }
    }
    else if( ! strcmp( fFieldName, "PixelShiftMod" ) || ! strcmp( fFieldName, "PixelShiftMod1" ) ) {  // JB,CB,PLR 2015/03/24
      // reading kShift for each device of this type
      for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
        pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod] = 3; // default value
        if( iMod>0 ) nextField();
        if( strstr( fFieldName, "PixelShiftMod" ) ) {
          read_item(pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod]);
        }
        else {
          printf( "WARNING in ReadDAQBoardParameters: field %s found while 'PixelShiftMod%d' expected!\n    -> Board %d of type %d lacks kShift.", fFieldName, iMod+1, iMod+1, aBoardNumber);
        }
      }
    }
    else
    {
      if (  strcmp( fFieldName, "StatisticCells") &&  strcmp( fFieldName, "CmsNoiseCut") &&  strcmp( fFieldName, "MaxNbOfHits") &&  strcmp( fFieldName, "Name") &&  strcmp( fFieldName, "AnalysisGoal") )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }      
    nextField();
    
  } while (  strcmp( fFieldName, "StatisticCells") &&  strcmp( fFieldName, "CmsNoiseCut") &&  strcmp( fFieldName, "MaxNbOfHits") &&  strcmp( fFieldName, "Name") &&  strcmp( fFieldName, "AnalysisGoal") );
  
  // Copy the parameters of the last specified input to the inputs 
  //  for which they were not specified.
  // This mechanism allows to declare N inputs while providing info just for one.
  // JB 2013/08/16
  for (Int_t aInp=inputCounter+1; aInp<pAcqModuleParameter[aBoardNumber].Inputs; aInp++) {
    printf(" updating input %d with input %d\n", aInp, inputCounter);
    pAcqModuleParameter[aBoardNumber].Channels[aInp] = pAcqModuleParameter[aBoardNumber].Channels[inputCounter];
    pAcqModuleParameter[aBoardNumber].Bits[aInp] = pAcqModuleParameter[aBoardNumber].Bits[inputCounter];
    pAcqModuleParameter[aBoardNumber].SigBits[aInp] = pAcqModuleParameter[aBoardNumber].SigBits[inputCounter];
    pAcqModuleParameter[aBoardNumber].NbOfFramesPerChannel[aInp] = pAcqModuleParameter[aBoardNumber].NbOfFramesPerChannel[inputCounter];
  }
  
  if( MimosaSimuSetupDebug) {
    cout << "  Using " << pAcqModuleParameter[aBoardNumber].Devices << " " << pAcqModuleParameter[aBoardNumber].Name << " with " << pAcqModuleParameter[aBoardNumber].Inputs << " inputs and " << pAcqModuleParameter[aBoardNumber].Bits << " bits(signif. " << pAcqModuleParameter[aBoardNumber].SigBits << ") and data files ";
    for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
      cout << " (" << iMod+1 << ")" << pAcqModuleParameter[aBoardNumber].DeviceDataFile[iMod];
    }
    cout << ", and PixelShift "; 
    for( Int_t iMod=0; iMod<pAcqModuleParameter[aBoardNumber].Devices; iMod++) {
      cout << " (" << iMod+1 << ")" << pAcqModuleParameter[aBoardNumber].PixelShiftMod[iMod];
    }
    cout << endl;
  }
    
  if ( ! strcmp( fFieldName, "Name") ) { // new board    
    ReadDAQBoardParameters( aBoardNumber+1 );
  } //end if new board
  

}
//==================================================================
void MimosaSimuSetup::ReadAnalysisParameters() 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter for Analysis 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified: JB 2013/06/21 add MaxTracksExGeom and ExGeomatrix params
  // Modified: JB 2013/07/17 new UserFlag param
  // Modified: JB 2013/09/21 new histoRange param
  // Modified: JB 2014/01/17 new AnalysisGoal param
  
  cout << endl << " - Reading Parameter for analysis" << endl;
  
  // Default values in case those parameters are set in the config file
  AnalysisParameter.CacheSize      = 0;
  AnalysisParameter.MaxNbOfHits    = 1000;
  AnalysisParameter.MinNbOfHits    = 0;
  AnalysisParameter.TrackChi2Limit = 1000.0;
  AnalysisParameter.Submatrices    = 0;
  for(int i=0;i<fMaxSubmatrices;i++) {
    AnalysisParameter.PixelSizeU[i]            = 0.0;
    AnalysisParameter.PixelSizeV[i]            = 0.0;
    AnalysisParameter.PixelsInRaw[i]           = 0;
    AnalysisParameter.PixelsInColumn[i]        = 0;
    AnalysisParameter.Matrixtype[i]            = -1;
    AnalysisParameter.MinSeedIndex[i]          = -1;
    AnalysisParameter.MaxSeedIndex[i]          = -1;
    AnalysisParameter.MinSeedCol[i]            = -1;
    AnalysisParameter.MaxSeedCol[i]            = -1;
    AnalysisParameter.MinSeedRow[i]            = -1;
    AnalysisParameter.MaxSeedRow[i]            = -1;
    AnalysisParameter.MaxNofPixelsInCluster[i] = 1000;
    AnalysisParameter.MinNofPixelsInCluster[i] = 1;
    AnalysisParameter.MinSeedCharge[i]         = -1.0e+6;
    AnalysisParameter.MinClusterCharge[i]      = -1.0e+6;
    AnalysisParameter.MinNeighbourCharge[i]    = -1.0e+6;
    AnalysisParameter.NoiseScope[i]            = 0.0;
    AnalysisParameter.Calibration[i]           = 1.0;
    AnalysisParameter.Geomatrices[i]           = 0;

    for(int j=0;j<fMaxGeomatrices;j++) {
      AnalysisParameter.Umin[i][j] = -1.0e+20;
      AnalysisParameter.Umax[i][j] = +1.0e+20;
      AnalysisParameter.Vmin[i][j] = -1.0e+20;
      AnalysisParameter.Vmax[i][j] = +1.0e+20;
    }
  }
  
  AnalysisParameter.MinHitsPerTrack = TrackerParameter.PlanesForTrackMinimum;
  AnalysisParameter.MaxTracksExGeom = -1;
  AnalysisParameter.ExGeomatrix = 0;
  AnalysisParameter.UserFlag = 0;
  AnalysisParameter.StatisticCells = 50;
  AnalysisParameter.CmsNoiseCut = 3.;
  AnalysisParameter.HistoChargeRange = 5000.; // e-, JB 2013/09/12
  AnalysisParameter.HistoSNRRange = 250;
  AnalysisParameter.HistoNoiseRange = 40; // e-
  AnalysisParameter.SavePlots = false;
  AnalysisParameter.DoTelResolutionMC = false;
  AnalysisParameter.MCEvents = 10000;
  AnalysisParameter.MCSeed   = 18383;
  AnalysisParameter.MCDoDisplay = false;
  AnalysisParameter.DoGaussianMS = true;
  AnalysisParameter.ResolutionScanSteps = 10;
  AnalysisParameter.ResolutionScanInit  =  2.0; // microns
  AnalysisParameter.ResolutionScanEnd   =  8.0; // microns
  AnalysisParameter.AnalysisGoal        = TString("");
  
  do {
    
    if ( ! strcmp( fFieldName, "StatisticCells") ) {
      read_item(AnalysisParameter.StatisticCells);
    }
    else if( ! strcmp( fFieldName, "CmsNoiseCut" ) ) {
      read_item(AnalysisParameter.CmsNoiseCut);
    }
    else if( ! strcmp( fFieldName, "MaxNbOfHits" ) ) {
      read_item(AnalysisParameter.MaxNbOfHits);
    }
    else if( ! strcmp( fFieldName, "MinNbOfHits" ) ) {
      read_item(AnalysisParameter.MinNbOfHits);
    }
    else if( ! strcmp( fFieldName, "TrackChi2Limit" ) ) {
      read_item(AnalysisParameter.TrackChi2Limit);
    }
    else if( ! strcmp( fFieldName, "SavePlots" ) ) {
      int temp;
      read_item(temp);
      if(temp == 1) AnalysisParameter.SavePlots = true;
      else          AnalysisParameter.SavePlots = false;
    }
    else if( ! strcmp( fFieldName, "DoTelResolutionMC" ) ) {
      int temp;
      read_item(temp);
      if(temp == 1) AnalysisParameter.DoTelResolutionMC = true;
      else          AnalysisParameter.DoTelResolutionMC = false;
      if(AnalysisParameter.DoTelResolutionMC) {
	cout << endl;
	cout << "Analysis Parameters:: Doing MC generation for evaluation of Telescope resolution at DUT!" << endl;
	cout << endl;
      }
    }
    else if( ! strcmp( fFieldName, "MCEvents" ) ) {
      read_item(AnalysisParameter.MCEvents);
    }
    else if( ! strcmp( fFieldName, "MCSeed" ) ) {
      read_item(AnalysisParameter.MCSeed);
    }
    else if( ! strcmp( fFieldName, "MCDoDisplay" ) ) {
      int temp;
      read_item(temp);
      if(temp == 1) AnalysisParameter.MCDoDisplay = true;
      else          AnalysisParameter.MCDoDisplay = false;
    }
    else if( ! strcmp( fFieldName, "DoGaussianMS" ) ) {
      int temp;
      read_item(temp);
      if(temp == 1) AnalysisParameter.DoGaussianMS = true;
      else          AnalysisParameter.DoGaussianMS = false;
    }
    else if( ! strcmp( fFieldName, "TrackChi2Limit" ) ) {
      read_item(AnalysisParameter.TrackChi2Limit);
    }
    else if( ! strcmp( fFieldName, "ResolutionScanSteps" ) ) {
      read_item(AnalysisParameter.ResolutionScanSteps);
    }
    else if( ! strcmp( fFieldName, "ResolutionScanInit" ) ) {
      read_item(AnalysisParameter.ResolutionScanInit);
    }
    else if( ! strcmp( fFieldName, "ResolutionScanEnd" ) ) {
      read_item(AnalysisParameter.ResolutionScanEnd);
    }
    else if( ! strcmp( fFieldName, "MinHitsPerTrack" ) ) {
      read_item(AnalysisParameter.MinHitsPerTrack);
    }
    else if( ! strcmp( fFieldName, "MaxTracksExGeom" ) ) {
      read_item(AnalysisParameter.MaxTracksExGeom);
    }
    else if( ! strcmp( fFieldName, "ExGeomatrix" ) ) {
      read_item(AnalysisParameter.ExGeomatrix);
    }
    else if( ! strcmp( fFieldName, "UserFlag" ) ) { // JB 2013/07/17
      read_item(AnalysisParameter.UserFlag);
    }
    else if( ! strcmp( fFieldName, "HistoChargeRange" ) ) { // JB 2013/09/12
      read_item(AnalysisParameter.HistoChargeRange);
    }
    else if( ! strcmp( fFieldName, "HistoSNRRange" ) ) { // JB 2013/09/12
      read_item(AnalysisParameter.HistoSNRRange);
    }
    else if( ! strcmp( fFieldName, "HistoNoiseRange" ) ) { // JB 2013/09/12
      read_item(AnalysisParameter.HistoNoiseRange);
    }
    else if( ! strcmp( fFieldName, "AnalysisGoal" ) ) { // JB 2014/01/17
      Char_t mmm[50];
      read_strings(mmm, 50);
      AnalysisParameter.AnalysisGoal = TString(mmm);
    }
    else if( ! strcmp( fFieldName, "Submatrices" ) || ! strcmp( fFieldName, "Subamtrices" )) { // note the second test to take into account typo propagated years ago !
      read_item(AnalysisParameter.Submatrices);
      if ( AnalysisParameter.Submatrices>fMaxSubmatrices ) {
        printf( "WARNING: required number of submatrices %d exceeds maximum allowed %d!\n  -> restricted to %d.\n", AnalysisParameter.Submatrices, fMaxSubmatrices, fMaxSubmatrices);
        AnalysisParameter.Submatrices = fMaxSubmatrices;
      }
    }
    else
    {
      if (  strcmp( fFieldName, "PixelSizeU" ) && !fConfigFileStream.eof() )
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }         
    nextField();
    
  } while (  strcmp( fFieldName, "PixelSizeU" ) && !fConfigFileStream.eof() );

  if( MimosaSimuSetupDebug) {
    cout << " - Final analysis cuts: " << endl;
    cout << "  Statistic cell size " << AnalysisParameter.StatisticCells << endl;
    cout << "  CMS cut " << AnalysisParameter.CmsNoiseCut << endl;
    cout << "  # hits max " << AnalysisParameter.MaxNbOfHits << endl;
    cout << "  # hits min " << AnalysisParameter.MinNbOfHits << endl;
    cout << "  Chi2 max " << AnalysisParameter.TrackChi2Limit << endl;
    if ( AnalysisParameter.MaxTracksExGeom > -1 ) cout << "  Max #tracks allowed in geomatrix " << AnalysisParameter.ExGeomatrix << " is " << AnalysisParameter.MaxTracksExGeom << endl;
    cout << "  # Submatrices " << AnalysisParameter.Submatrices << endl;
    cout << endl;
  }

  
}
//==================================================================
void MimosaSimuSetup::ReadSubmatrixParameters( Int_t aSubmatrixNumber) 
{
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter for Submatrix 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //
  // JB 2013/01/16
  // Modified JB 2013/07/17 Matrixtype introduced
  // Modified: JB 2013/08/21,22 new SeedIndex/Col/Row limits
  // Modified: JB 2013/09/12 new params MinNofPixelsInCluster
  // Modified: JB 2013/11/08 new cuts on charges
  // Modified: JB 2014/01/21 new cut on cluster charge

  if ( aSubmatrixNumber<0 || fMaxSubmatrices<=aSubmatrixNumber ) {
    printf( "WARNING: trying to add submatrix %d beyond the maximum number allowed (%d)\n  --> nothing done!\n", aSubmatrixNumber+1, fMaxSubmatrices);
    return;
  }
  else {
    cout << "  * adding submatrix " << aSubmatrixNumber << endl;
  }
  
  // some initialization
  AnalysisParameter.Matrixtype[aSubmatrixNumber] = 1;
  AnalysisParameter.Calibration[aSubmatrixNumber] = 1.;
  AnalysisParameter.Geomatrices[aSubmatrixNumber] = 0;
  AnalysisParameter.MinSeedIndex[aSubmatrixNumber] = 0; // JB 2013/08/21
  AnalysisParameter.MaxSeedIndex[aSubmatrixNumber] = 0;
  AnalysisParameter.MinSeedCol[aSubmatrixNumber] = 0; // JB 2013/08/22
  AnalysisParameter.MaxSeedCol[aSubmatrixNumber] = 0;
  AnalysisParameter.MinSeedRow[aSubmatrixNumber] = 0; // JB 2013/08/22
  AnalysisParameter.MaxSeedRow[aSubmatrixNumber] = 0;
  AnalysisParameter.MinNofPixelsInCluster[aSubmatrixNumber] = 1; // JB 2013/09/12
  AnalysisParameter.MinSeedCharge[aSubmatrixNumber] = -1000.; // JB 2013/11/08
  AnalysisParameter.MinClusterCharge[aSubmatrixNumber] = -1000.; // JB 2014/01/21
  AnalysisParameter.MinNeighbourCharge[aSubmatrixNumber] = -1000.; // JB 2013/11/08

  do {
    if( ! strcmp( fFieldName, "PixelSizeU" ) ) {
      read_item(AnalysisParameter.PixelSizeU[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "PixelSizeV" ) ) {
      read_item(AnalysisParameter.PixelSizeV[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "PixelsInRaw" ) ) {
      read_item(AnalysisParameter.PixelsInRaw[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "PixelsInColumn" ) ) {
      read_item(AnalysisParameter.PixelsInColumn[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinSeedIndex" ) ) { // JB 2013/08/21
      read_item(AnalysisParameter.MinSeedIndex[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MaxSeedIndex" ) ) { // JB 2013/08/21
      read_item(AnalysisParameter.MaxSeedIndex[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinSeedCol" ) ) { // JB 2013/08/21
      read_item(AnalysisParameter.MinSeedCol[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MaxSeedCol" ) ) { // JB 2013/08/21
      read_item(AnalysisParameter.MaxSeedCol[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinSeedRow" ) ) { // JB 2013/08/22
      read_item(AnalysisParameter.MinSeedRow[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MaxSeedRow" ) ) { // JB 2013/08/22
      read_item(AnalysisParameter.MaxSeedRow[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MaxNofPixelsInCluster" ) ) {
      read_item(AnalysisParameter.MaxNofPixelsInCluster[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinNofPixelsInCluster" ) ) { // JB 2013/09/12
      read_item(AnalysisParameter.MinNofPixelsInCluster[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinSeedCharge" ) ) { // JB 2013/11/08
      read_item(AnalysisParameter.MinSeedCharge[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinClusterCharge" ) ) { // JB 2014/01/21
      read_item(AnalysisParameter.MinClusterCharge[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "MinNeighbourCharge" ) ) { // JB 2013/11/08
      read_item(AnalysisParameter.MinNeighbourCharge[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "Matrixtype" ) ) { // JB 2013/07/17
      read_item(AnalysisParameter.Matrixtype[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "Calibration" ) ) {
      read_item(AnalysisParameter.Calibration[aSubmatrixNumber]);
    }
    else if( ! strcmp( fFieldName, "NoiseScope" ) ) {
      read_item(AnalysisParameter.NoiseScope[aSubmatrixNumber]); // JB 2010/09/06
    }
    else if( strstr( fFieldName, "GeoMatrix" ) || strstr( fFieldName, "Geomatrix" ) ) {
      if( AnalysisParameter.Geomatrices[aSubmatrixNumber] == fMaxGeomatrices ) {
        printf( "WARNING in ReadSubmatrixParameters: maximum nb of Geomatrix reached, taking no more!\n");
      }
      else {
        read_item(AnalysisParameter.Umin[aSubmatrixNumber][AnalysisParameter.Geomatrices[aSubmatrixNumber]]);
        nextItem(':');
        read_item(AnalysisParameter.Umax[aSubmatrixNumber][AnalysisParameter.Geomatrices[aSubmatrixNumber]]);
        nextItem(':');
        read_item(AnalysisParameter.Vmin[aSubmatrixNumber][AnalysisParameter.Geomatrices[aSubmatrixNumber]]);
        nextItem(':');
        read_item(AnalysisParameter.Vmax[aSubmatrixNumber][AnalysisParameter.Geomatrices[aSubmatrixNumber]]);
        AnalysisParameter.Geomatrices[aSubmatrixNumber]++;
      }
    }
    else
    {
      if (  strcmp( fFieldName, "PixelSizeU" ))
      {
        cout << "WARNING : parameter '" << fFieldName << "' in config file is not understood !" << endl;
        getRidOfLine();
      }
    }
  
    nextField();
  
  } while (  strcmp( fFieldName, "PixelSizeU" ) && !fConfigFileStream.eof() );

  if( MimosaSimuSetupDebug )  {
    cout << "  * Submatrix " << aSubmatrixNumber << endl;
    cout << "    pixel size U " << AnalysisParameter.PixelSizeU[aSubmatrixNumber] << " V " << AnalysisParameter.PixelSizeV[aSubmatrixNumber] << endl;
    cout << "    matrix type " << AnalysisParameter.Matrixtype[aSubmatrixNumber] << endl;
    cout << "    # pixels " << AnalysisParameter.PixelsInRaw[aSubmatrixNumber] << " in raw, " << AnalysisParameter.PixelsInColumn[aSubmatrixNumber] << " in column" << endl;
    cout << "    min - max # pixels in hit " << AnalysisParameter.MinNofPixelsInCluster[aSubmatrixNumber] << " - " << AnalysisParameter.MaxNofPixelsInCluster[aSubmatrixNumber] << endl;
    cout << "    range of seed index: " << AnalysisParameter.MinSeedIndex[aSubmatrixNumber] << " to " << AnalysisParameter.MaxSeedIndex[aSubmatrixNumber] << endl; 
    cout << "    calibration " << AnalysisParameter.Calibration[aSubmatrixNumber] << endl;
    cout << "    noiseScope " << AnalysisParameter.NoiseScope[aSubmatrixNumber] << endl;
    cout << "    # Geomatrices " << AnalysisParameter.Geomatrices[aSubmatrixNumber] << endl;
      for( Int_t il=0; il<AnalysisParameter.Geomatrices[aSubmatrixNumber]; il++) {
        cout << "    GeoMatrix " << il << " " << AnalysisParameter.Umin[aSubmatrixNumber][il] << " < U < " << AnalysisParameter.Umax[aSubmatrixNumber][il] << "; " << AnalysisParameter.Vmin[aSubmatrixNumber][il] << " < V < " << AnalysisParameter.Vmax[aSubmatrixNumber][il] << endl;
      }
    }

  if ( ! strcmp( fFieldName, "PixelSizeU" ) ) {
    ReadSubmatrixParameters( aSubmatrixNumber+1);
  }
  
}
//==================================================================
void MimosaSimuSetup::ReadConfiguration() 
{
  // Reads telescope plane and detector configuration from configuration file
  // The call order of each method matters!
  // 1) Run param.
  // 2) Tracker param.
  // 2bis) Geometrical Parameters for the IVI vertexing
  // 3) Plane and/or Ladder param.
  // 4) DAQ param.
  // 5) DAQ Module/Board param.
  // 6) Final analysis param.
  //
  // Modified: JB 2013/01/16, use new parsing methods
  // Modified: VR 2016/06/30, run number is get from DSession and DataSubDirPrefixXXX (XXX = run #) concatenated to DataPath if set
  
  cout << endl << " -*-*- MimosaSimuSetup User Constructor -*-*- " << endl;  

  // -+-+-+-+-+--+-+-
  // --- Initialization
  // -+-+-+-+-+--+-+-
  printf(" - Reading Setup from %s\n", fConfigPathAndFileName.Data());
  
  fAddedPlanes = 0;
  fAddedLadders = 0;
  
  // -+-+-+-+-+--+-+-
  // --- open config file:
  // -+-+-+-+-+--+-+-
  fConfigFileStream.open(fConfigPathAndFileName);
  Bool_t answer = fConfigFileStream.fail();
  if (answer && gROOT->IsBatch()) {
    printf("ERROR ! Can't read file %s\nQuit (because in batch mode)\n", fConfigPathAndFileName.Data());
    gSystem->Exit(-1);
  }
  while(answer) {
    cout << "enter correct file name \n";
    cin >> fConfigPathAndFileName;
    printf(" - Reading Setup from %s\n", fConfigPathAndFileName.Data());
    fConfigFileStream.open(fConfigPathAndFileName);
    answer=fConfigFileStream.fail();
  } 
  //TString command = TString("cat ") + fConfigPathAndFileName;
  //gSystem->Exec(command.Data());
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Run Parameter 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  ReadRunParameters();
  //RunParameter.Number = fSession->GetRunNumber();// VR 20014/06/30 Set the run number from DSession(fRunNumber), which is set by MimosaAnalysis::InitSession(const Int_t TheRun)
  //if ( strcmp(RunParameter.DataSubDirPrefix,""))// VR 20014/06/30 
  //if DataSubDirPrefix is given, concatenate it with the run number to the DataPath
  //{
    //sprintf( RunParameter.DataPath, "%s/%s%d/", RunParameter.DataPath,RunParameter.DataSubDirPrefix, RunParameter.Number ); // VR 2014/06/30
  //}
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Tracker 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  ReadTrackerParameters();

  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Experiment Geometry
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  //ReadExperimentGeometryParameters(); //Valerian 2015/02/02@18h50 : if this method causes bug, tell me so I can correct bug (In a previous commit it was comment ?!)
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Planes 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  pPlaneParameter  = new PlaneParameter_t[TrackerParameter.Planes];
  ChannelUse      = new Int_t*[TrackerParameter.Planes];
  if ( ! strcmp( fFieldName, "Inputs" ) ) ReadPlaneParameters(fAddedPlanes); // contains iterative call to subsequent planes
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameters of the Ladders 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  pLadderParameter  = new LadderParameter_t[TrackerParameter.Ladders];
  ChannelUse      = new Int_t*[TrackerParameter.Ladders];
  if ( ! strcmp( fFieldName, "LadderID" ) ) ReadLadderParameters(fAddedLadders); // contains iterative call to subsequent ladders
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter of the Data Acquisition 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  ReadDAQParameters();
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter of the Data Acquisition boards in this run
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  pAcqModuleParameter = new AcqModuleParameter_t[AcqParameter.ModuleTypes];
  ReadDAQBoardParameters(0); // contains iterative call to subsequent Module
  
  
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  // Parameter for Analysis 
  // -+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-
  ReadAnalysisParameters();
  if ( ! strcmp( fFieldName, "PixelSizeU" ) ) {
    ReadSubmatrixParameters( 0); // contains iterative call to subsequent submatrices
  }
  else {
    printf( "WARNING in ReadAnalysisParameter: Missing submatrix parameters (PixelSize, ...)!\n");  
  }

  cout << endl << " -*-*- MimosaSimuSetup User Constructor DONE -*-*- " << endl;

  return;
  
}
//================================================================== 
MimosaSimuSetup::MimosaSimuSetup(const MimosaSimuSetup& c)
{
  copy(c);
}
//==================================================================  
MimosaSimuSetup& MimosaSimuSetup::operator=(const MimosaSimuSetup& a)
{
  copy(a);
  return *this;
}
//==================================================================
void MimosaSimuSetup::copy(const MimosaSimuSetup& a)
{
  // Prepare a Copy into this class
}
//================================================================== 
void MimosaSimuSetup::nextItem(Char_t delimiter)
{
  // Move the file parsing pointer to the character "delimiter"
  
  Char_t c;
  do {
    fConfigFileStream >> c;
    if (MimosaSimuSetupDebug>1)  cout << c;
  } while (c != delimiter);
}
//================================================================== 
void MimosaSimuSetup::nextField()
{
  // Move the file parsing pointer to the character "delimiter"
  // store in fFieldName the previous chain of characters 
  //  (skipping spaces, new lines and the delimiter)
  //
  // Modified JB 2012/12/20 include recording of fieldName
  // Modified JB 2013/01/16 skip comment lines
  
  Char_t delimiter = ':';
  Int_t k = 0;
  Char_t c, previousC;
  do {
    fConfigFileStream >> c;
    //cout << "|" << c << "|";
    if( c != '\n' && c != ' ' && c != delimiter && c != '.') {
      fFieldName[k] = c;
      previousC = c;
    }
    if( (c == '\n') || (c == ' ') || (k==fFieldMaxLength-1) || (c == '.'))
    {
      k = 0;
    }
    else if (c != delimiter) {
      k++;
    }
    if ( c=='/' && previousC=='/' ) {
      getRidOfLine();
      k = 0;
    }
    //cout << "k=" << k ;
  } while (c != delimiter && !fConfigFileStream.eof() );
  fFieldName[k]='\0';
  if (MimosaSimuSetupDebug>1)  cout << "field = " << fFieldName << endl;
}
//==================================================================
void MimosaSimuSetup::read_r3(G4ThreeVector &arg)
{
  // Modified BH 2013/08/21 memory leak removed
  
  Double_t co[3] = {0., 0., 0.}; // BH 2013/08/21
  for (Int_t k = 0; k < 3; k++) {
    if( k>0 ) nextItem(':'); // already positionned for 1st value
    fConfigFileStream >> co[k];
    if (MimosaSimuSetupDebug>1) cout << co[k] << endl;
    arg.setX(co[0]);
    arg.setY(co[1]);
    arg.setZ(co[2]);
  }
}
//==================================================================
void MimosaSimuSetup::read_item(Int_t &arg)
{
  //nextItem(':');
  fConfigFileStream >> arg;
  if (MimosaSimuSetupDebug>1){
    cout << "value = " << arg << endl;
  }
}
//==================================================================
void MimosaSimuSetup::read_item(UInt_t &arg)
{
  //nextItem(':');
  fConfigFileStream >> arg;
  if (MimosaSimuSetupDebug>1){
    printf("value =%d/%x\n",arg,arg); 
   }
}
//==================================================================
void MimosaSimuSetup::read_item(Float_t &arg)
{

  // reads values from configuration file

  //nextItem(':');
  fConfigFileStream >> arg;
  if (MimosaSimuSetupDebug>1) cout << "value = " << arg << endl;
}
//==================================================================
void MimosaSimuSetup::read_strings(Char_t *aString, Int_t aLength)
{

  // reads a string of a given max length from configuration file
  // The strings is expected to be contained by double quotes : "
  // JB 2009/05/25

  Int_t k = 0;
  Char_t c;
  // First, go to the " delimiter
  //nextItem('"');
  do {
    fConfigFileStream >> c;
    //cout << c;
  } while (c != '"');  
  // Now, read the value up to the next " delimiter
  do {
    fConfigFileStream >> c;
    //cout << c;
    if ((c != '"') && (k < aLength)) {
      aString[k] = c;
      k++;
    }
  } while (c != '"');
  aString[k]='\0'; // end properly the string when it is shorter than max length
  if (MimosaSimuSetupDebug>1) cout << "value = " << aString << endl;
  
}
//==================================================================
void MimosaSimuSetup::read_TStrings(TString& TheString, Int_t aLength)
{

  Char_t aString[500];
  // reads a string of a given max length from configuration file
  // The strings is expected to be contained by double quotes : "
  // JB 2009/05/25

  Int_t k = 0;
  Char_t c;
  // First, go to the " delimiter
  //nextItem('"');
  do {
    fConfigFileStream >> c;
    //cout << c;
  } while (c != '"');  
  // Now, read the value up to the next " delimiter
  do {
    fConfigFileStream >> c;
    //cout << c;
    if ((c != '"') && (k < aLength)) {
      aString[k] = c;
      k++;
    }
  } while (c != '"');
  aString[k]='\0'; // end properly the string when it is shorter than max length
  if (MimosaSimuSetupDebug>1) cout << "value = " << aString << endl;

  TheString = TString(aString);
  
}
//==================================================================
void MimosaSimuSetup::getRidOfLine()
{
  
  // Simply get rid of all character till the line ends,
  // line is not expected to exceeds 200 charaters.
  // 
  // JB, 2013/01/16
  
  Char_t line[250];
  fConfigFileStream.getline( line, 250);

}
//==================================================================
void  MimosaSimuSetup::CheckExpSetupParameters(void)
{
  

  //Test of the experimental setup parameters
  cout << endl;
  if(TrackerParameter.ExpSetup == TString("") || 
     (TrackerParameter.ExpSetup == TString("Beam-Test") || TrackerParameter.ExpSetup == TString("beam-test"))
    ) {
    CheckBeamTestParameters();
  }
  else if(CheckIfSourceExists(TrackerParameter.ExpSetup)) {
    CheckSourceParameters();
  }
  else {
    cout << endl;
    cout << "The experimental setuo parameter ExpSetup needs to be specified to one of the following values." << endl;
    cout << " - Beam-Test,"       << endl;
    for(int isource=0;isource<int(_Source_list.size());isource++) cout << " - " << _Source_list[isource].Data() << "," << endl;
    cout << "Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }
  
  return;
  
}
//==================================================================
void  MimosaSimuSetup::CheckBeamTestParameters(void)
{
  
  if(TrackerParameter.ExpSetup == TString("")) {
     TrackerParameter.ExpSetup = TString("Beam-Test");
     cout << "ExpSetup  parameter not specified. Setting it to default value:  " << TrackerParameter.ExpSetup.Data() << endl;
  }
  else if(TrackerParameter.ExpSetup == TString("Beam-Test") || TrackerParameter.ExpSetup == TString("beam-test")) {
    TrackerParameter.ExpSetup = TString("Beam-Test");
    cout << "ExpSetup  parameter set to user value:  " << TrackerParameter.ExpSetup.Data() << endl;
  }
    
  if(TrackerParameter.BeamType == TString("")) {
    TrackerParameter.BeamType = TString("pion-");
    cout << "    BeamType             parameter not specified. Setting it to defaul value:  " << TrackerParameter.BeamType.Data() << endl;
  }
  else cout << "    BeamType             parameter set to user value: " << TrackerParameter.BeamType.Data() << endl;
    
    
  if(TrackerParameter.BeamNparticles == -1.0) {
    TrackerParameter.BeamNparticles = 1.0;
    cout << "    BeamNparticles       parameter not specified. Setting it to defaul value:  " << TrackerParameter.BeamNparticles << endl;
  }
  else if(TrackerParameter.BeamNparticles < 0.0) {
    TrackerParameter.BeamNparticles = 1.0;
    cout << "    BeamNparticles       parameter set to negative value. Setting it to defaul value:  " << TrackerParameter.BeamNparticles << endl;
  }
  else cout << "    BeamNparticles       parameter set to user value: " << TrackerParameter.BeamNparticles << endl;
    
  if(TrackerParameter.BeamRandNparticles) {
    cout << "        Generate number of particle with Poisson distribution with mean value :  " << TrackerParameter.BeamNparticles << endl;
  }
  else cout << "        Generate exactly " << TrackerParameter.BeamNparticles << " particles per event" << endl;
    
  if(TrackerParameter.BeamMomentum == -1.0) {
    TrackerParameter.BeamMomentum = 120.0;
    cout << "    BeamMomentum         parameter not specified. Setting it to defaul value: " << TrackerParameter.BeamMomentum << " GeV/c" << endl;
  }
  else if(TrackerParameter.BeamMomentum < 0) {
    TrackerParameter.BeamMomentum = 120.0;
    cout << "    BeamMomentum         parameter set to user value is negative. Setting it to defaul value: " << TrackerParameter.BeamMomentum << " GeV/c" << endl;
  }
  else cout << "    BeamMomentum         parameter set to user value: " << TrackerParameter.BeamMomentum << " GeV/c" << endl;
    
  if(TrackerParameter.BeamDirection == G4ThreeVector(0.0,0.0,1.0)) {
    cout << "    BeamDirection        parameter not specified. Setting it to defaul value: (" << TrackerParameter.BeamDirection(0) << "," << TrackerParameter.BeamDirection(1) << "," << TrackerParameter.BeamDirection(2) << ")" << endl;
  }
  else {
    //Normalizing Beam direction vector to unit module
    if(TMath::Abs(TrackerParameter.BeamDirection.mag()) < 1.0e-8) {
      G4cerr << "BeamDirection    vector has zero module. Check your inputs. Exiting now!!!" << G4endl;
      assert(false);
    }
    TrackerParameter.BeamDirection /= TrackerParameter.BeamDirection.mag();
    cout << "    BeamDirection        parameter set to user value: (" << TrackerParameter.BeamDirection(0) << "," << TrackerParameter.BeamDirection(1) << "," << TrackerParameter.BeamDirection(2) << ")" << endl;
  }
    
  if(TrackerParameter.BeamAngularSpread == G4ThreeVector(0.0,0.0,1.0)) {
    cout << "    BeamAngularSpread    parameter not specified. Setting it to defaul value: (" << TrackerParameter.BeamAngularSpread(0) << "," << TrackerParameter.BeamAngularSpread(1) << ")" << endl;
  }
  else if(TrackerParameter.BeamAngularSpread(0) < 0.0 || TrackerParameter.BeamAngularSpread(1) < 0.0) {
    TrackerParameter.BeamAngularSpread.setX(0.0);
    TrackerParameter.BeamAngularSpread.setY(0.0);
    cout << "    BeamAngularSpread    parameter has negative value. Setting it to defaul value: (" << TrackerParameter.BeamAngularSpread(0) << "," << TrackerParameter.BeamAngularSpread(1) << ")" << endl;
  } 
  else cout << "    BeamAngularSpread    parameter set to user value: (" << TrackerParameter.BeamAngularSpread(0) << "," << TrackerParameter.BeamAngularSpread(1) << ")" << endl;
    
  if(TrackerParameter.BeamMomentumSpread == -1) {
    TrackerParameter.BeamMomentumSpread = 1.0e-10;
    cout << "    BeamMomentumSpread   parameter not specified. Setting it to defaul value: " << TrackerParameter.BeamMomentumSpread << endl;
  }
  else if(TrackerParameter.BeamMomentumSpread < 0.0) {
    TrackerParameter.BeamMomentumSpread = 0.0;
    cout << "    BeamMomentumSpread   parameter has negative value. Setting it to defaul value: " << TrackerParameter.BeamMomentumSpread << endl;
  }
  else cout << "    BeamMomentumSpread   parameter set to user value: " << TrackerParameter.BeamMomentumSpread << endl;
    
  if(TrackerParameter.BeamOrigin == G4ThreeVector(0.0,0.0,-1.0e+6) || TrackerParameter.BeamOriginSpread(0) < 0.0 || TrackerParameter.BeamOriginSpread(1) < 0.0) {
    cout << "    BeamOrigin or BeamOriginSpread  parameters not specified. Will perform uniform generation in lowest Z face of the world volume cube" << endl;
  }
  else {
    cout << "    BeamOrigin and BeamOriginSpread  parameters specified by user. Will perform Gaussian generation at BeamOrigin position with BeamOriginSpread." << endl;
    cout << "    BeamOrigin           parameter set to user value: (" <<  TrackerParameter.BeamOrigin(0)/1000 << "," << TrackerParameter.BeamOrigin(1)/1000 << "," << TrackerParameter.BeamOrigin(2)/1000 << ") mm" << endl;
    cout << "    BeamOriginSpread     parameter set to user value: (" <<  TrackerParameter.BeamOriginSpread(0)/1000 << "," << TrackerParameter.BeamOriginSpread(1)/1000 << ") mm" << endl;
  }
  
  return;
  
}
//==================================================================
void  MimosaSimuSetup::CheckSourceParameters(void)
{
  
  if(TrackerParameter.SourcePosition == G4ThreeVector(-999.0,-999.0,-999.0)) {
    cout << "    SourcePosition         parameter not specified. Exiting now!!!" << endl;
    assert(false);
  }
  else cout << "    SourcePosition         parameter set to user value: (" << TrackerParameter.SourcePosition(0) << "," << TrackerParameter.SourcePosition(1) << "," << TrackerParameter.SourcePosition(2) << ") mm" << endl;
  
  if(TrackerParameter.SourceTilt == G4ThreeVector(-999.0,-999.0,-999.0)) {
    TrackerParameter.SourceTilt == G4ThreeVector(0.0,0.0,0.0);
    cout << "    SourceTilt             parameter not specified. Setting it to default value: (" << TrackerParameter.SourceTilt(0) << "," << TrackerParameter.SourceTilt(1) << "," << TrackerParameter.SourceTilt(2) << ") deg" << endl;
  }
  else cout << "    SourceTilt             parameter set to user value: (" << TrackerParameter.SourceTilt(0) << "," << TrackerParameter.SourceTilt(1) << "," << TrackerParameter.SourceTilt(2) << ") deg" << endl;
  
  if(TrackerParameter.SourceRadius == -1.0) {
    cout << "    SourceRadius           parameter not specified. Exiting now!!!" << endl;
    assert(false);
  }
  else cout << "    SourceRadius           parameter set to user value: " << TrackerParameter.SourceRadius << " mm" << endl;
  
  if(TrackerParameter.SourceActivity == -1.0) {
    cout << "    SourceActivity         parameter not specified. Exiting now!!!" << endl;
    assert(false);
  }
  else cout << "    SourceActivity         parameter set to user value: " << TrackerParameter.SourceActivity << " decays/sec" << endl;
  
  if(TrackerParameter.SourceSensorROTime == -1.0) {
    cout << "    SourceSensorROTime     parameter not specified. Exiting now!!!" << endl;
    assert(false);
  }
  else cout << "    SourceSensorROTime     parameter set to user value: " << TrackerParameter.SourceSensorROTime << " usec" << endl;
  
  return;
  
}
//==================================================================
void  MimosaSimuSetup::InitSourceList(void)
{

  _Source_list.clear();
  _Source_list.push_back("Source-Sr90");
  _Source_list.push_back("Source-Fe55");
  
  return;
  
}
//==================================================================
bool  MimosaSimuSetup::CheckIfSourceExists(TString SourceName)
{
  
  bool ItExists = false;
  for(int isource=0;isource<int(_Source_list.size());isource++) {
    if(_Source_list[isource] == SourceName) {
      ItExists = true;
      break;
    }
  }
  
  return ItExists;
  
}
//==================================================================
void  MimosaSimuSetup::CheckPlaneDigitizationParameters(Int_t aPlaneNumber)
{
  
  if(pPlaneParameter[aPlaneNumber].PlaneDigitization == TString("")) {
    cout << " WARNING: NO DIGITIZATION WILL BE PERFORMED for plane " << aPlaneNumber+1 << endl;
    return;
  }
  
  if(pPlaneParameter[aPlaneNumber].AnalysisMode < 2) {
    //Strip readout
    cout << "  WARNING: Not able to perform digitization for analysis mode < 2. No digitization will be performed!!!" << endl;
    pPlaneParameter[aPlaneNumber].PlaneDigitization = TString("");
    
    return;
    
  }
  else if(pPlaneParameter[aPlaneNumber].AnalysisMode == 2) {
    //Pixel readout with analog output
    cout << "Performing digitization with analogue output for plane " << aPlaneNumber+1 << ". Using model " << pPlaneParameter[aPlaneNumber].PlaneDigitization << endl;    
    
    if((pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise   > 0.0)     && 
       (pPlaneParameter[aPlaneNumber].PlaneDigitizeCalib   > 0.0)     && 
       (pPlaneParameter[aPlaneNumber].PlaneDigitizeADCbits > 0)       && 
       (pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin  != -999.0) &&
       (pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax  != -999.0) && 
       (pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin  <  pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax)) {
      cout << "  Using the following digitization parameters:" << endl;
      cout << "    - Sensor noise (electrons)        = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise   << endl;
      cout << "    - Sensor calibration (elec/Volts) = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeCalib   << endl;
      cout << "    - Sensor number of ADC bits       = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCbits << endl;
      cout << "    - Sensor ADC dynamic range        = (" << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin << "," << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax << ") volts" << endl;
    }
    else {
      cout << "  WARNING: The below parameters either not set of set to wrong alues:" << endl;
      cout << "    - Sensor noise (PlaneDigitizeNoise)                                   = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise << " (electrons)" << endl;
      cout << "    - Sensor calibration (PlaneDigitizeCalib)                             = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeCalib << " (elec/Volt)" << endl;
      cout << "    - Sensor number of ADC bits (PlaneDigitizeADCbits)                    = "  << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCbits << endl;
      cout << "    - Sensor ADC dynamic range  (PlaneDigitizeADCMin,PlaneDigitizeADCMax) = (" << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMin << "," << pPlaneParameter[aPlaneNumber].PlaneDigitizeADCMax << ") volts" << endl;
      cout << "  NO DIGITIZATION WILL BE PERFORMED!!!" << endl;
      pPlaneParameter[aPlaneNumber].PlaneDigitization = TString("");
    }
    
    return;
    
  }
  else if(pPlaneParameter[aPlaneNumber].AnalysisMode == 3) {
    //Pixel readout with digital output
    cout << "Performing digitization with digital output for plane " << aPlaneNumber+1 << ". Using model " << pPlaneParameter[aPlaneNumber].PlaneDigitization << endl;
    
    if(pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise > 0.0 && 
       pPlaneParameter[aPlaneNumber].PlaneDigitizeOcc   > 0.0 && 
       pPlaneParameter[aPlaneNumber].PlaneDigitizeThre  > 0.0) {
      cout << "  Using the following digitization parameters:" << endl;
      cout << "    - Sensor noise (electrons)       = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise << endl;
      cout << "    - Threshold (x noise)            = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeThre  << endl;
      cout << "    - Sensor average Noise occupancy = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeOcc   << endl;
    }
    else {      
      cout << "  WARNING: The below parameters either not set of set to wrong alues:" << endl;
      cout << "    - Sensor noise (electrons)       = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeNoise << endl;
      cout << "    - Threshold (x noise)            = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeThre  << endl;
      cout << "    - Sensor average Noise occupancy = " << pPlaneParameter[aPlaneNumber].PlaneDigitizeOcc   << endl;
      cout << "  NO DIGITIZATION WILL BE PERFORMED!!!" << endl;
      pPlaneParameter[aPlaneNumber].PlaneDigitization = TString("");
    }
    
  }
  
  return;
  
}
//==================================================================
  