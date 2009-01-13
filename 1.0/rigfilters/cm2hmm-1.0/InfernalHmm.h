/*
This file copyright (c) 2004, Zasha Weinberg
All rights reserved.

Redistribution and use in source and binary forms, 
with or without modification, are permitted 
provided that the following conditions are met:

- Redistributions of source code must retain the 
above copyright notice, this list of conditions 
and the following disclaimer. 
- Redistributions in binary form must reproduce 
the above copyright notice, this list of 
conditions and the following disclaimer in the 
documentation and/or other materials provided 
with the distribution. 
- Neither the name of the University of Washington 
nor the names of its contributors may be used to 
endorse or promote products derived from this 
software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// InfernalHmm class implements an HMM based on the representation of a CovarianceModel

enum Cm2Hmm_HmmBuildType {
	HmmBuildType_Original=0,  // also called "Compacted" or "type 0"
	HmmBuildType_separateMPandMLMR=1, // also called "Expanded" or "type 1".  for MATP node, each HMM node (left and right) has 2 ML states, 1 D state, 1 IL state
	HmmBuildType_separateMPMLMRD=2 // Didn't work well.  for MATP node, each HMM node (left and right) has 2 ML states, 2 D states, 1 IL state
};
extern void VerifyValidHmmBuildType (Cm2Hmm_HmmBuildType hmmType);
extern int GetNumExtraStatesVersusOriginalHmmBuildType (Cm2Hmm_HmmBuildType hmmType);
extern char *GetHmmBuildTypeDescription (Cm2Hmm_HmmBuildType hmmType);
extern Cm2Hmm_HmmBuildType GetHmmBuildTypeByText (const char *text);

enum HmmFileFormat {
	HmmFileFormat_Binary, // DEPRECATED: binary format - bad, mainly since it's not compatible between little-endian and big-endian CPUs, and also probably not between 32-bit and 64-bit CPUs
	HmmFileFormat_Text, // newer textual format (which is what the 'SaveInBinary' function saves by default)
	HmmFileFormat_HMMER // HMMER format.  Only works for compacted-type profile HMMs (expanded-type HMMs aren't technically profile HMMs, so they won't work with HMMER)
};

typedef VectorByCmState<std::list<CovarianceModel::State> > Hmm2CmStateVector;

// each CM state corresponds to 1 or 2 HMM states
struct Cm2HmmState {
	CovarianceModel::State hmmLeftState,hmmRightState; // -1 means not applicable
};
typedef VectorByCmState<Cm2HmmState> Cm2HmmStateVector;


class InfernalHmm : public CovarianceModelBase {
private:
	// these functions are not meaningful for HMMs (or at least not in the same way), so I'd like to prevent myself from using them by making them private
	float GetEndsc (State state) const;
	float GetBeginsc (State state) const;
	void SetBeginsc (State state,float sc);

protected:
	// to handle the endsc part of local alignments, we need this link info.
	// index endscLinksToLeftVector by an hmm state.  If this state is a left hmm state, then you'll get an empty list.  If it's a right hmm state, you'll get all the left hmm states that come from the same CM state as the right state.  For handling endsc in the hmm, we're allowed to skip from right to left at cost endsc.
	enum {MAX_ENDSC_LINK_TO_LEFTS=2};
	struct LinkToLeftInfo {
		State hmmLeftState;
		float endsc;
	};
	typedef FixedArrayWithSize<LinkToLeftInfo,MAX_ENDSC_LINK_TO_LEFTS> EndscLinksToLeft;
	typedef VectorByCmState<EndscLinksToLeft> EndscLinksToLeftVector;
	EndscLinksToLeftVector endscLinksToLeftVector;

	struct OtherStateInfo {
		bool isLeftState;
		bool isRightState;

		// we can have beginsc that correspond to left/right HMM nodes.  Sometimes a left/right HMM node can be related to 2 CM nodes, in which case e.g. it's left score will be the max of the 2 scores.  I'm worried about HMM nodes that are both left & right nodes, so to make things easier, I'm just going to keep the leftward & rightward scores totally separate.
		float leftwardBeginsc,rightwardBeginsc;
	};
	typedef VectorByCmState<OtherStateInfo> OtherStateInfoVector;
	OtherStateInfoVector otherStateInfoVector;
	Hmm2CmStateVector hmm2CmStateVector;
	Cm2HmmStateVector cm2HmmStateVector;

	struct ParentAndMyChildNum {
		State parentState;
		int myChildNum;
	};
	typedef std::list<ParentAndMyChildNum> ParentAndMyChildNumVector;
	struct NonSavedInfo {
		ParentAndMyChildNumVector parentAndMyChildNumVector;
		int childNumOfSelfLoop; // of -1 if this state has no self loop
	};
	typedef VectorByCmState<NonSavedInfo> NonSavedInfoVector;
	NonSavedInfoVector nonSavedInfoVector;

	std::string fromCmFileName;
	std::string fullBuildDescription;
	static int GetCurrFileFormatNum (void);
	int loadedFileFormat;

	Cm2Hmm_HmmBuildType hmmBuildType; // how were we built (redundant with fullBuildDescription, but nice to have to computer know for algorithms that require a particular method)

public:

	InfernalHmm (void);

	void AddBuildDescription (const std::string addedBuildDescription);
	void SetHmmBuildType (Cm2Hmm_HmmBuildType hmmBuildType_);
	Cm2Hmm_HmmBuildType GetHmmBuildType (void) const;
	void SetFromCmFileName (const std::string& _fromCmFileName);
	std::string GetBuildDescription (void) const;
	void SaveInDeprecatedBinaryFormat (const char *fileName);
	const std::string& GetFromCmFileName (void) const;
	int GetFileFormat (void) const;

	// by default, we don't make nodes for the HMMs, but if we want to save them in Infernal format, we need to make them up.  Also changes PASSTHRU_st --> S_st, which hopefully Infernal will be okay with
	void ConvertToInfernalSavableFormat (const CovarianceModel& cm);

	void HmmNormalize (void);

	// returns -1 if given state has no self-loop
	// caller must call BuildNonSavedInfoIfNecessary at some point before calling this, or have called CopyFrom on an HMM that had it called on (or was CopyFrom'd one that... etc)
	int GetChildNumOfSelfLoop (State state) const;

	void BuildNonSavedInfoIfNecessary (void);

	void Init (int numStates);

	void DumpInfernalHmm (FILE *file,const CovarianceModel& cm) const;
	bool LoadInBinary (FILE *cmFile);
	void SaveInBinary (FILE *cmFile);
	bool LoadInBinary (const char *cmFileName);
	void SaveInBinary (const char *cmFileName);
	void SaveInDeprecatedBinaryFormat (FILE *file); // really do it
	void CopyReverseOf(const InfernalHmm& t);
	void CopyFrom (const InfernalHmm& t);
	void CopyFromWithEscHack(const InfernalHmm& t);
	void MirrorFromWithHackedExtraInfo(CM_t *cm);
	void ClobberIR (void); // for creating HMMs from cmbuild; kill the IR in the start state
	void SaveInFormat (const char *hmmFileName,HmmFileFormat hmmFileFormat);
	static HmmFileFormat GetDefaultHmmFileFormat (void);

	void SetHmm2CmState (const Hmm2CmStateVector& _hmm2CmStateVector);
	void SetCm2HmmState (const Cm2HmmStateVector& _cm2HmmStateVector);

	// hackquery
	void AllocHmmData (void);

	void SetLeftState (State state,bool isLeftState);
	void SetRightState (State state,bool isRightState);
	void SetLeftwardBeginsc (State state,float sc);
	void SetRightwardBeginsc (State state,float sc);

	inline int GetNumEndscLinksToLeft (State state) const {
		return endscLinksToLeftVector[state].size();
	}
	inline State GetEndscLinkToLeft_State (State state,int linkNum) const {
		return endscLinksToLeftVector[state][linkNum].hmmLeftState;
	}
	inline float GetEndscLinkToLeft_Endsc (State state,int linkNum) const {
		return endscLinksToLeftVector[state][linkNum].endsc;
	}

	inline void SetNumEndscLinksToLeft (State state,int n) {
		endscLinksToLeftVector[state].resize(n);
	}
	inline void SetEndscLinkToLeft_State (State state,int linkNum,State hmmLeftState) {
		endscLinksToLeftVector[state][linkNum].hmmLeftState=hmmLeftState;
	}
	inline void SetEndscLinkToLeft_Endsc (State state,int linkNum,float endsc) {
		endscLinksToLeftVector[state][linkNum].endsc=endsc;
	}

	inline bool IsLeftState (State state) const {
		return otherStateInfoVector[state].isLeftState;
	}
	inline bool IsRightState (State state) const {
		return otherStateInfoVector[state].isRightState;
	}
	inline float GetLeftwardBeginsc (State state) const {
		return otherStateInfoVector[state].leftwardBeginsc;
	}
	inline float GetRightwardBeginsc (State state) const {
		return otherStateInfoVector[state].rightwardBeginsc;
	}

	inline const CovarianceModel::StateList& GetCmState (State hmmState) const {
		return hmm2CmStateVector[hmmState];
	}
	inline const Hmm2CmStateVector& GetCmStateVector (void) const {
		return hmm2CmStateVector;
	}
	CovarianceModel::NodeList GetCmNode (const CovarianceModel& cm,State hmmState) const;

	inline State GetHmmLeftStateOfCmState (CovarianceModel::State cmState) const {
		return cm2HmmStateVector[cmState].hmmLeftState;
	}
	inline State GetHmmRightStateOfCmState (CovarianceModel::State cmState) const {
		return cm2HmmStateVector[cmState].hmmRightState;
	}
	inline const Cm2HmmStateVector GetHmmStateVector (void) const {
		return cm2HmmStateVector;
	}

	// note: this will not necessarily return states that are unique to a given node, just states that implement CM nodes.  There will in general be states that relate to more than one nodes at splice points from splicing bifurcations together.
	int GetNumStatesForCmNode (const CovarianceModel& cm,CovarianceModel::Node cmNode) const;

	struct StateFromNode {
		State state;
		bool isRight; // when coming from a Node, a state is only either left or right
		bool isNormal; // as opposed to post-insert state
	};
	typedef std::list<StateFromNode> StateFromNodeList;
	void GetHmmStatesByCmNode (StateFromNodeList& get_hmmStates,const CovarianceModel& cm,CovarianceModel::Node cmNode) const;
	void GetNormalHmmStatesOfLeftOrRightByCmNode (StateList& get_stateList,const CovarianceModel& cm,CovarianceModel::Node cmNode,bool getIsRight) const;
	struct EdgeInfo {
		State fromState,toState;
		int childNum; // relative to fromState, of course
		bool isRightSide; // else on left side of HMM
	};
	typedef std::list<EdgeInfo> EdgeInfoList;
	void GetHmmEdgesByCmNode (EdgeInfoList& get_edgeInfoList,const StateFromNodeList& hmmStates,const CovarianceModel& cm,CovarianceModel::Node cmNode);

	void DumpProbabilitiesCsvByNode(FILE *nodeDump,const CovarianceModel& cm,bool doValues); // doValues as opposed to do headers for csv
	// copies them into this Hmm
	void CopyProbabilitiesForNode(const InfernalHmm& sourceHmm,const CovarianceModel& cm,const CovarianceModel::Node cmNode)
;

	bool AreSingletEmissionScoresSameForAllNucs(const State& state) const;

	void DumpCmAndHmmCorrespondingScores (const char *outFileName,const CovarianceModel& cm);

	typedef float PairInflationMatrix[MAXABET][MAXABET];
	void ComputePairInflationMatrix (PairInflationMatrix& pairInflationMatrix,const CovarianceModel& cm,CovarianceModel::State cmState) const;

	typedef float SingletEmissionScores[MAXABET];
	// starts with scores for the given infernalHmm (this pointer).  adjusts them to take into account only those pairs that are PENALIZED in boolBasePairMatrix
	void AdjustLeftAndRightEscForRewardedPairsOnly (SingletEmissionScores& leftNucScores,SingletEmissionScores& rightNucScores,const CovarianceModel& cm,const CovarianceModel::State state,const BoolBasePairMatrix& boolBasePairMatrix) const;
};
