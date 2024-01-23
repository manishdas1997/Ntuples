/**
 *  @copyright Copyright 2019 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.h
 */

#ifndef EVENTCATEGORIZER_H
#define EVENTCATEGORIZER_H

#include <JPetCommonTools/JPetCommonTools.h>
#include <JPetUserTask/JPetUserTask.h>
#include "../ModularDetectorAnalysis/EventCategorizerTools.h"
#include "reconstructor.h"
#include <JPetEvent/JPetEvent.h>
#include <Hits/JPetBaseHit/JPetBaseHit.h>
#include <Hits/JPetMCRecoHit/JPetMCRecoHit.h>
#include <Hits/JPetPhysRecoHit/JPetPhysRecoHit.h>
//#include <JPetRawMCHit/JPetRawMCHit.h>
#include <JPetTimer/JPetTimer.h>
#include <vector>
#include <map>
#include <TCutG.h>

using namespace std;
class JPetWriter;

#ifdef __CINT__
#define override
#endif

static const double kLightVelocity_cm_ns = 29.9792458;

enum EvtType{

  UNKNOWN = 1,
  OPS,
  OPS_AND_1SCATTER,
  B2B_AND_SCATTER,
  B2B_AND_PROMPT,
  OPS_AND_PROMPT,
  OPS_AND_SCATTER_PROMPT,
  OPS_SCATTER_AND_PROMPT,
  B2B_AND_2SCATTER,
  OPS_AND_2SCATTER,
  B2B_SCATTER_AND_PROMPT,
  B2B_PROMPT_AND_SCATTER,
  PROMPT_SCATTER,
  B2B_PRMT_2SCAT,
  PRMT_B2B_2SCAT,
  ALL_THREE_FROM_SAME_VTX,
  TWO_FROM_SAME_VTX,
  NONE_FROM_SAME_VTX,
  OTHER,
  OPS_SCATTER_NONSENSITIVE,
  B2B_SCATTER_NONSENSITIVE,
  PMT_SCATTER_NONSENSITIVE
  
};
typedef std::pair< double, double > coordinates;


class EventCategorizer : public JPetUserTask{
public:
	EventCategorizer(const char * name);
	virtual ~EventCategorizer(){};
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;	

	double meanTime_3Hits;
	
	double Magnitude(TVector3 p1);
	void saveEvents(const std::vector<JPetEvent>& event);
	
	void initialiseHistograms();

	//vector<const JPetMCRecoHit*> tot( const JPetEvent &event);
	//vector<const JPetMCRecoHit*> scatTest( vector<const JPetMCRecoHit*> hits_tot );

	vector<const JPetPhysRecoHit*> tot( const JPetEvent *event);
	vector<const JPetPhysRecoHit*> scatTest( vector< const JPetPhysRecoHit*> hits_tot );
	
	//double OPS_events( const JPetEvent *event, const JPetTimeWindowMC& time_window_mc );


	//void genOPSCalculation( vector<JPetBaseHit*> hits,  const JPetTimeWindowMC& time_window_mc, EvtType& event_type);                                                                    
        //double genHitsCalculation( vector<JPetBaseHit*> hits, const JPetMCRecoHit &mchit1, const JPetMCRecoHit& mchit2, const JPetMCRecoHit& mchit3, EvtType& event_type, const JPetTimeWindow MC& time_window_mc);
	//<JPetBaseHit*> hits, const JPetMCRecoHit &mchit1, const JPetMCRecoHit& mchit2, const JPetMCRecoHit& mchit3, EvtType& event_type, const JPetTimeWindowMC& time_window_mc);

	
        void ReconHitsCalculation( vector<const JPetPhysRecoHit*> hits, EvtType& event_type,const JPetTimeWindowMC& time_window_mc);

	TVector3 calculateAnnihilationPoint(  const JPetPhysRecoHit* hit1,   const JPetPhysRecoHit* hit2);
	static double calculateTOF(  const JPetPhysRecoHit* hitA,   const JPetPhysRecoHit* hitB);
	static double calculateTOF(double time1, double time2);

	
	double scatterTest( vector<const JPetPhysRecoHit*> hits, EvtType event_type );
     	
	std::vector<double> TOT(std::vector<double> eng);
       
	
	double calculateAngle(const JPetPhysRecoHit* hit1,  const JPetPhysRecoHit* hit2 );
	double calculatedLOR( vector<const JPetPhysRecoHit*> hits, EvtType event_type);
	
	std::vector<double> calcLOR_recsAnnhPt( vector<const JPetPhysRecoHit*> hits, TVector3 decayPt);
	
	//double calculateAngleSum (JPetBaseHit* hit1, JPetBaseHit* hit2, JPetBaseHit* hit3);
	
	vector<double> distance_hits (vector<const JPetPhysRecoHit*> hits);
	std::pair<double, double> hitDistance_tDiff (vector<const JPetPhysRecoHit*> hits);
	
	std::vector<double> hitsDis_diff ( vector<const JPetPhysRecoHit*> hits);
	double calculateAngle_3D( const JPetPhysRecoHit* hit1, const JPetPhysRecoHit* hit2 );	
	
	//void check3gSameVtx( vector<JPetBaseHit*> hits, const JPetTimeWindowMC& time_window_mc);
	
	TVector3 annhPoint ( const JPetPhysRecoHit* hit1,  const JPetPhysRecoHit* hit2, const JPetTimeWindowMC& time_window_mc );
	double angularRes_rec (   const JPetPhysRecoHit* hit1,   const JPetPhysRecoHit* hit2, const JPetTimeWindowMC& time_window_mc );
	//double angularRes_gen ( const JPetBaseHit* hit1, const JPetBaseHit* hit2, const JPetTimeWindowMC& time_window_mc );
	

	TVector3 cal2gAnnhPoint( const JPetEvent *event, EvtType event_type );
	bool ellipse_cut_values( double angle_sum, double dlor_min);

	double decayTime_3hits(std::vector<const JPetPhysRecoHit*> hits , EvtType& event_type);
	//double lifeTime(const JPetEvent *event , EvtType& event_type, std::vector<double> prompt_emTime, int size);

	void lifeTime( std::vector<double> promptTime, std::vector<std::pair<double,EvtType>> hitTime );

	
	double operator_ops_gen, operator_ops_recs;
	double operator_opsScat_gen, operator_opsScat_recs;
        

private:
	//double decayTime, pmtEmmTime;
	//std::vector<double> prompt_time, ops_time;
	
   
protected:

        Reconstructor * fReconstructor;

	const std::string kMC = "Save_MC_bool";
	bool fIsMC = true;

	std::map<EvtType, std::string> EventTypeNames_MC= {
	  {UNKNOWN, "3 hit evts"},
	  {OPS, "oPs"},
	  {OPS_AND_1SCATTER, "oPs and scattered photon"},
	  {B2B_AND_SCATTER, "back-to-back and scattered photon"},
	  {B2B_AND_PROMPT, "back-to-back and prompt photon"},
	  {OPS_AND_PROMPT, "oPs and prompt photon"},
	  {OPS_AND_SCATTER_PROMPT, "o-Ps, prompt, sactter prompt"},
	  {OPS_SCATTER_AND_PROMPT, "o-Ps, scatter, prompt"},
	  {OPS_AND_2SCATTER, "oPs and 2 scattered photon"},
	  {B2B_AND_2SCATTER, "back-to-back and 2 scattered photon"},
	  {B2B_SCATTER_AND_PROMPT, "b2b, b2b scattered and prompt"},
	  {B2B_PROMPT_AND_SCATTER, "b2b, prompt and prompt scatter"},
	  {PROMPT_SCATTER, "prompt and 2 scattered prompt"},
	  {PRMT_B2B_2SCAT, "prmt, b2b, b2b-doubleScat"},
	  {B2B_PRMT_2SCAT, "b2b, prmt, prmt-doubleScat"},
	  {ALL_THREE_FROM_SAME_VTX, "3g from same vtx"},
	  {TWO_FROM_SAME_VTX, "Any two from same vtx"},
	  {NONE_FROM_SAME_VTX, "all from different vtx"},
	  {OTHER, "other evts"},
	  {OPS_SCATTER_NONSENSITIVE, "ops scatter nSensitive"},
	  {B2B_SCATTER_NONSENSITIVE, "b2b scatter nSensitive"},
	  {PMT_SCATTER_NONSENSITIVE, "pmt scatter nSensitive"}
	  
	};

	std::map<EvtType, std::string> EventTypeNames_data = {
	  {UNKNOWN, "3 hit evts"}
	};
	
	EvtType event_type;
	std::map<EvtType, std::string> EventTypeNames;

	double pmtEmmTime;
   
	std::vector<double> prompt_time;

	std::vector<std::pair<double, EvtType>> threeHitTime;
	
       	
	  
};

#endif /* !EVENTCATEGORIZER_H */
